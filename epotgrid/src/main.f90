! /***************************************
!   epotgrid 
!  22 Aug. 2017  written by D.Kawata
! ****************************************/

program enez
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
      use grid_particle
      use mext
   
      implicit none
      include 'mpif.h'

      integer id,i,j,k,pn,flag,ns,ng,np,istep,ngpt,ngp
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer nstep,step,is,ic
      integer nval
      integer flagcom,flagr,flagmext,flaggp,flagrot,flagrlog
      integer ierr
! for z=0 grid data
      integer nr,nth,nz,ir,ith,iz,iths,ithe,ip,nthi
      integer nx,ny,ix,iy,iys,iye,nyi
      double precision rran(0:1),zran(0:1)
      double precision lnri,lnro,dlnr,dz,dth,dr,ri,ro
      double precision thpi,zpi,rpi,potri
      double precision xran(0:1),yran(0:1)
      double precision xpi,ypi,dx,dy
      double precision hmp
      double precision fr,z0
! rotation degree
      double precision rotang
! particle info
      double precision hpm,r2dp,thp,zplane
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! ***  for filename ***
      character sname*6,filename*60,stepname*6,fileo*60
      integer logs,slogs
      double precision fs
      character soname*3
      integer uni,count
! for work
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:),tx(:),ty(:)

! ***** MPI Initialization *****
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      if(myrank.eq.0) then
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif
      if(nprocs.gt.NCPU) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: increase NCPU=',NCPU
        endif
        call MPI_FINALIZE()
        stop
      endif

      call setkernel()

! *****   Open Input File   ******
      if(myrank.eq.0) then
        open(50,file='./ini/input.dat',status='old')
        read(50,*) nstep,step
        read(50,*) flaggp,flagrlog
        if(flaggp.eq.0) then
           read(50,*) nr,nth,nz
           read(50,*) rran(0),rran(1)
           read(50,*) zran(0),zran(1)
        else
           read(50,*) nx,ny,nz
           read(50,*) xran(0),xran(1)
           read(50,*) yran(0),yran(1)
           read(50,*) zran(0),zran(1)
        endif
        read(50,*) flagmext
        read(50,*) hmp
        read(50,*) flagrot
        read(50,*) rotang
        close(50)
        if(nstep.gt.1) then
          write(6,*) ' nstep=',nstep
        else
          write(6,*) ' step=',step 
        endif
        if(flaggp.eq.0) then
          write(6,*) ' radial circular grid'
          write(6,*) ' grid at z=0 nR,nth,nz=',nr,nth,nz
          write(6,*) ' R range (kpc)=',rran(0),rran(1)
          if(flagrlog.ne.0) then
            write(6,*) ' R grid in log scale'
          else
            write(6,*) ' R grid in linear scale'
          endif
        else
          write(6,*) ' square grid'
          write(6,*) ' grid at z=0 nx,ny,nz=',nx,ny,nz
          write(6,*) ' x range (kpc)=',xran(0),xran(1)
          write(6,*) ' y range (kpc)=',yran(0),yran(1)
        endif
        write(6,*) ' z range (kpc)=',zran(0),zran(1)
! kpc -> 100 kpc (GCD+ unit)
        if(flaggp.eq.0) then
          rran(0)=rran(0)/LUKPC
          rran(1)=rran(1)/LUKPC
        else
          xran(0)=xran(0)/LUKPC
          xran(1)=xran(1)/LUKPC
          yran(0)=yran(0)/LUKPC
          yran(1)=yran(1)/LUKPC
        endif
        zran(0)=zran(0)/LUKPC
        zran(1)=zran(1)/LUKPC
        write(6,*) ' smoothing length (kpc)=',hmp
        hmp=hmp/LUKPC
        if(flagmext.ne.0) then
          write(6,*) ' add external potential from ini/mext.dat'
        endif
        if(flagrot.ne.0) then
          write(6,*) ' rotate anti-clockwise (deg) =',rotang
        endif
      endif

! send the data to the other node
      call MPI_BCAST(flaggp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(flagrlog,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      nval=7
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        if(flaggp.eq.0) then
          tivr(2)=nr
          tivr(3)=nth
        else
          tivr(2)=nx
          tivr(3)=ny
        endif
        tivr(4)=nz
        tivr(5)=flagmext
        tivr(6)=flagrot
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      if(flaggp.eq.0) then
        nr=tivr(2)
        nth=tivr(3)
      else
        nx=tivr(2)
        ny=tivr(3)
      endif
      nz=tivr(4)
      flagmext=tivr(5)
      flagrot=tivr(6)
      if((nprocs.gt.nth.and.flaggp.eq.0).or.(nprocs.gt.ny.and.flaggp.ne.0)) then
        if(myrank.eq.0) then
          write(6,*) ' Error: nprocs > nth'
        endif
        call MPI_ABORT()
      endif

      deallocate(tivr)

      if(flaggp.eq.0) then
        nval=6
      else
        nval=8
      endif
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=hmp
        tdvr(1)=rotang
        if(flaggp.eq.0) then
          tdvr(2)=rran(0)
          tdvr(3)=rran(1)
          tdvr(4)=zran(0)
          tdvr(5)=zran(1)
        else
          tdvr(2)=xran(0)
          tdvr(3)=xran(1)
          tdvr(4)=yran(0)
          tdvr(5)=yran(1)
          tdvr(6)=zran(0)
          tdvr(7)=zran(1)
        endif
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      hmp=tdvr(0)      
      rotang=tdvr(1)      
      if(flaggp.eq.0) then
        rran(0)=tdvr(2)
        rran(1)=tdvr(3)
        zran(0)=tdvr(4)
        zran(1)=tdvr(5)
      else
        xran(0)=tdvr(2)
        xran(1)=tdvr(3)
        yran(0)=tdvr(4)
        yran(1)=tdvr(5)
        zran(0)=tdvr(6)
        zran(1)=tdvr(7)
      endif
      deallocate(tdvr)

! set how to read the data
      flagr=1

      if(flagmext.ne.0) then
! set up the external potential
        call setmext()
      endif

      if(nstep.gt.1.and.myrank.eq.0) then
        open(49,file='./ini/file.dat',status='old')
      endif
      do istep=1,nstep
        if(nstep.gt.1.and.myrank.eq.0) then
          read(49,*) step
        endif    
        call MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! *** read the data ***
        call rdata(step,ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns,ai,tu,flagr)
        npt=ngt+nst
        np=ng+ns
        if(myrank.eq.0) then
           write(6,*) ' step, t(GYR)=',step,tu*TMUGYR
           write(6,*) ' ngt,ndmt,ndm1t,nst=',ngt,ndmt,ndm1t,nst
        endif

! rotation
        if(flagrot.ne.0) then
! deg to radian
          rotang=rotang*M_PI/180.0d0
! for baryon
          if(np.gt.0) then
            allocate(tx(0:np-1))
            allocate(ty(0:np-1))
            do i=0,np-1
              tx(i)=x_p(i)*dcos(-rotang)-y_p(i)*dsin(-rotang)
              ty(i)=x_p(i)*dsin(-rotang)+y_p(i)*dcos(-rotang)
              x_p(i)=tx(i)
              y_p(i)=ty(i)
            enddo
            deallocate(tx)
            deallocate(ty)
! for test output
            open(60,file='rots.dat',status='unknown')
            do i=ng,np-1,10
              pn=list_ap(i)
              write(60,'(3(1pE13.5))') x_p(pn),y_p(pn),z_p(pn)
            enddo
          endif
! for DM
          if(ndm.gt.0) then
            allocate(tx(0:ndm-1))
            allocate(ty(0:ndm-1))
            do i=0,ndm-1
              tx(i)=x_dm(i)*dcos(-rotang)-y_dm(i)*dsin(-rotang)
              ty(i)=x_dm(i)*dsin(-rotang)+y_dm(i)*dcos(-rotang)
              x_dm(i)=tx(i)
              y_dm(i)=ty(i)
            enddo
            deallocate(tx)
            deallocate(ty)
          endif
        endif
! domain decomposition
        if(ngt.gt.0.or.nst.gt.0) then
          call ddecb(npt,np,ng,ns)
! Note: ng and ns are not set, and no longer true
        endif
        if(ndmt.gt.0) then
          call ddecdm(ndmt,ndm)

!          write(fileo,'(a3,i3.3)') 'dmp',myrank
!          open(60,file=fileo,status='unknown')
!          do i=0,ndm-1
!            write(60,'(4(1pE13.5))') x_dm(i),y_dm(i),z_dm(i),h_dm(i)
!          enddo
!          close(60)

        endif

        if(myrank.eq.0) then
          write(6,*) ' after domain dec np,ndm=',np,ndm
        endif

! generate grid particle
        if(flaggp.eq.0) then
          ngpt=nr*nth*nz
          if(flagrlog.eq.0) then
! linear R grid
            lnri=0.0d0
            lnro=0.0d0
            dlnr=0.0d0
            ri=rran(0)
            ro=rran(1)
            if(nr.gt.1) then
              dr=(ro-ri)/dble(nr-1)
            else
              dr=0.0d0
            endif
          else
            ri=0.0d0
            ro=0.0d0
            dr=0.0d0
            lnri=dlog(rran(0))
            lnro=dlog(rran(1))
            if(nr.gt.2) then
              dlnr=(lnro-lnri)/dble(nr-2)
            else
              dlnr=0.0d0
            endif
          endif
          if(nz.gt.1) then
            dz=(zran(1)-zran(0))/dble(nz-1)
          else
            dz=0.0d0 
          endif
          dth=2.0d0*M_PI/dble(nth)
          if(myrank.eq.0) then
            write(6,*) ' dz=z range /(Nz-1), grid stars from zran0'
            write(6,*) ' dth=2 Pi/Nth, grid stars from 0 deg'
            if(flagrlog.eq.0) then
              write(6,*) ' dr=r range/(nr-1), grid starts from rran0'
              write(6,*) ' dR, dth, dz=',dr,dth,dz
            else
              write(6,*) &
                ' dlnR=ln R range /(Nr-2), grid starts from 0 and rran0'
              write(6,*) ' dlnR, dth, dz=',dlnr,dth,dz
            endif
          endif

          call para_range(0,nth-1,nprocs,myrank,iths,ithe)

          nthi=ithe-iths+1
          ngp=nr*nz*nthi
          if(myrank.eq.0) then
            write(6,*) ' number of grids total and myrank 0 =',ngpt,ngp
            write(6,*) ' ith range for myrank 0 =',iths,ithe
          endif

          call allocate_gridparticle(ngp,nr,nz,nthi)

          z0=dabs(zran(1))
          iz0=0
          ip=0
          do ith=iths,ithe
            thpi=dth*dble(ith)
            do iz=0,nz-1
              zpi=zran(0)+dz*dble(iz)
! to find z=0
              if(dabs(zpi).lt.z0) then
                iz0=iz
                z0=dabs(zpi)
              endif
              do ir=0,nr-1
! store the particle position
                pngrid(ir,iz,ith-iths)=ip
                if(flagrlog.eq.0) then
                  rpi=ri+dble(ir)*dr
                else
                  if(ir.eq.0) then
! ir=0 grid point at r=1
                    rpi=0.0d0
                  else
! ir=1 grid point at r=rran(0)
                    rpi=dexp(lnri+dlnr*dble(ir-1))
                  endif
                endif
                xp(ip)=rpi*dcos(thpi)
                yp(ip)=rpi*dsin(thpi)
                zp(ip)=zpi
                hp(ip)=hmp
                epotp(ip)=0.0d0
                dvxp(ip)=0.0d0
                dvyp(ip)=0.0d0
                dvzp(ip)=0.0d0
                ip=ip+1
              enddo
            enddo
          enddo
          if(ip.ne.ngp) then
            write(6,*) ' Error: ip not equal ngp at myrank, ip, ngp=',myrank &
             ,ip,ngp
          endif
          if(myrank.eq.0) then
            write(6,*) ' z=0 grid iz=',iz0,z0
          endif
        else
          ngpt=nx*ny*nz
          if(nx.gt.1) then
            dx=(xran(1)-xran(0))/dble(nx-1)
          else
            dx=0.0d0
          endif
          if(ny.gt.1) then
            dy=(yran(1)-yran(0))/dble(ny-1)
          else
            dy=0.0d0
          endif
          if(nz.gt.1) then
            dz=(zran(1)-zran(0))/dble(nz-1)
          else
            dz=0.0d0 
          endif
          dth=2.0d0*M_PI/dble(nth)
          if(myrank.eq.0) then
            write(6,*) ' dxyz=xyzz range /(Nzxyz-1), grid stars from xyzran'
            write(6,*) ' dx, dy, dz=',dx,dy,dz
          endif

! parallel alog y-axis.
          call para_range(0,ny-1,nprocs,myrank,iys,iye)

          nyi=iye-iys+1
          ngp=nx*nz*nyi
          if(myrank.eq.0) then
            write(6,*) ' number of grids total and myrank 0 =',ngpt,ngp
            write(6,*) ' iy range for myrank 0 =',iys,iye
          endif

          call allocate_gridparticle(ngp,nx,nz,nyi)

          z0=dabs(zran(1))
          iz0=0
          ip=0
          do iy=iys,iye
            ypi=yran(0)+dy*dble(iy)
            do iz=0,nz-1
              zpi=zran(0)+dz*dble(iz)
! to find z=0
              if(dabs(zpi).lt.z0) then
                iz0=iz
                z0=dabs(zpi)
              endif
              do ix=0,nx-1
! store the particle position
                pngrid(ix,iz,iy-iys)=ip
                xpi=xran(0)+dx*dble(ix)
                xp(ip)=xpi
                yp(ip)=ypi
                zp(ip)=zpi
                hp(ip)=hmp
                epotp(ip)=0.0d0
                dvxp(ip)=0.0d0
                dvyp(ip)=0.0d0
                dvzp(ip)=0.0d0
                ip=ip+1
              enddo
            enddo
          enddo
          if(ip.ne.ngp) then
            write(6,*) ' Error: ip not equal ngp at myrank, ip, ngp=',myrank &
             ,ip,ngp
          endif
          if(myrank.eq.0) then
            write(6,*) ' z=0 grid iz=',iz0,z0
          endif
        endif
 
        if(npt.gt.0) then
          call treebuild(np)
        endif
        if(ndmt.gt.0) then
          call dmtreebuild(ndm)
        endif

        call treeforce(ngp,npt,ndmt)

        if(flagmext.ne.0) then
          do i=0,ngp-1
            rpi=dsqrt(xp(i)**2+yp(i)**2+zp(i)**3)
            if(rpi.gt.10.0d0**SI_lro) then
              potri=mextr(SI_nmext-1)/rpi
              epotp(i)=epotp(i)-G*potri
            else if(rpi.gt.0.0d0) then
              ir=int((dlog10(rpi)-SI_lri)/SI_dlr)
              if(ir.lt.0) then
                ir=0
              else if(ir.ge.SI_nmext-1) then
                ir=SI_nmext-2
              endif
              potri=(mextr(ir)+(rpi-rmext(ir)) &
                *(mextr(ir+1)-mextr(ir))/(rmext(ir+1)-rmext(ir)))/rpi
              epotp(i)=epotp(i)-potri
            endif
            if(rpi.gt.0.0d0) then
              potri=potri/(rpi**2)
            else
              potri=0.0d0
            endif
            dvxp(i)=dvxp(i)-xp(i)*G*potri
            dvyp(i)=dvyp(i)-yp(i)*G*potri
            dvzp(i)=dvzp(i)-zp(i)*G*potri
          enddo
        endif

!        write(fileo,'(a6,i3.3)') 'gridpf',myrank
!        open(60,file=fileo,status='unknown')
!        do i=0,ngp-1
!          rpi=dsqrt(xp(i)**2+yp(i)**2)
!          fr=(dvxp(i)*xp(i)+dvyp(i)*yp(i))/rpi
!          write(60,'(9(1pE13.5))') xp(i),yp(i),zp(i),epotp(i) &
!            ,dvxp(i),dvyp(i),dvzp(i),dsqrt(-fr*rpi)*VUKMS,rpi*LUKPC
!        enddo
!        close(60)

        if(flaggp.eq.0) then
          call output(step,ngp,flaggp,flagrlog &
            ,nr,nz,nth,nthi,rran,zran,0.0,0.0)
        else
          call output(step,ngp,flaggp,flagrlog &
            ,nx,nz,ny,nyi,xran,zran,yran)
        endif

      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program 
