! /***************************************
!   pslice ver.1
!  14 Aug. 2017  written by D.Kawata
! ****************************************/

program pslice
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
      use smap_particle
   
      implicit none
      include 'mpif.h'

      integer id,i,j,k,pn,flag,ns,ng,np,istep,nsp,flagr
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer nstep,step,is,ic
      integer nval
      integer flag3d,flagrot,flagsp,flagid,ids,ide
      integer ierr
! particle info
      double precision,allocatable :: temp_p(:),nh_p(:)
! grid info
      integer nxgrid,nygrid,nsmapval
      double precision zsl
      double precision xrange(0:1),yrange(0:1),zrange(0:1)
      double precision hmp,hmpx,hmpy,pdx,pdy,xpm,ypm,xlsm,xhsm,ylsm,yhsm
      double precision,allocatable :: fval(:,:,:)
      double precision tmass
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! ***  for filename ***
      character sname*6,filename*60,stepname*6,fileo*60
! for work
      integer,allocatable :: slist(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:) &
        ,tx(:),ty(:),tz(:),tm(:),rp(:)

! ***** MPI Initialization *****
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      if(myrank.eq.0) then
        write(6,*) ' pslice ver. 1 14/08/17'
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif
      if(nprocs.gt.NCPU) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: increase NCPU=',NCPU
        endif
        call MPI_FINALIZE()
        stop
      endif

! number of values
      nsmapval=8

! *****   Open Input File   ******
      if(myrank.eq.0) then
        open(50,file='./ini/input.dat',status='old')
        read(50,*) nstep,step
        read(50,*) nxgrid,nygrid
        read(50,*) flag3d
        read(50,*) xrange(0),xrange(1)
        read(50,*) yrange(0),yrange(1)
        read(50,*) zrange(0),zrange(1)
        read(50,*) zsl
        read(50,*) flagrot
        read(50,*) flagsp
        read(50,*) hmp
        read(50,*) flagid
        read(50,*) ids,ide
        close(50)
        close(50)
        if(nstep.gt.1) then
          write(6,*) ' nstep=',nstep
        else
          write(6,*) ' step=',step 
        endif
        if(flagsp.eq.0) then
          write(6,*) ' analysis for gas'
!        else if(flagsp.eq.1) then
!          write(6,*) ' analysis for DM, nbm=',nbm
!        else if(flagsp.eq.2) then      
!          write(6,*) ' analysis for stars, nbm=',nbm
        else 
          write(6,*) ' flagsp should be 0, i.e. only gas.',flagsp
          stop
        endif
        if(flag3d.eq.0) then
          write(6,*) ' grid nx, ny=',nxgrid,nygrid
          write(6,*) ' x range (kpc)=',xrange(0),xrange(1)
          write(6,*) ' y range (kpc)=',yrange(0),yrange(1)
          write(6,*) ' at z=',zsl
        else if(flag3d.lt.0) then
          write(6,*) 'nx,nz = ',nxgrid,nygrid
          write(6,*) ' x range = ',xrange(0),xrange(1)
          write(6,*) ' z range = ',yrange(0),yrange(1)
          write(6,*) ' at y = ',zsl
        else 
          zsl=0.0
          write(6,*) ' 2D case'
          write(6,*) ' selected zrange=',zrange(0),zrange(1)
        endif
        if(flag3d.le.0) then
          zrange(0)=0.0
          zrange(1)=0.0
        endif
        if(flagrot.ne.0) then
          write(6,*) ' rotation axis set to -z'
        endif
        write(6,*) ' smoothing length if gas, minimum (kpc)=',hmp
        if(flagid.ne.0) then
          write(6,*) ' only in id range =',ids,ide
        endif
! send the data to the other node
      endif

      nval=9
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        tivr(2)=nxgrid
        tivr(3)=nygrid
        tivr(4)=flag3d
        tivr(5)=flagrot
        tivr(6)=flagsp
        tivr(7)=flagid
        tivr(8)=ids
        tivr(9)=ide
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      nxgrid=tivr(2)
      nygrid=tivr(3)
      flag3d=tivr(4)
      flagrot=tivr(5)
      flagsp=tivr(6)
      flagid=tivr(7)
      ids=tivr(8)
      ide=tivr(9)

      deallocate(tivr)

      nval=8
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=xrange(0)
        tdvr(1)=xrange(1)
        tdvr(2)=yrange(0)
        tdvr(3)=yrange(1)
        tdvr(4)=zrange(0)
        tdvr(5)=zrange(1)
        tdvr(6)=hmp
        tdvr(7)=zsl
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      xrange(0)=tdvr(0)
      xrange(1)=tdvr(1)
      yrange(0)=tdvr(2)
      yrange(1)=tdvr(3)
      zrange(0)=tdvr(4)
      zrange(1)=tdvr(5)
      hmp=tdvr(6)
      zsl=tdvr(7)

      deallocate(tdvr)

! set kernel
      call setkernel(flag3d)

! smoothing length for x and y
      pdx=((xrange(1)-xrange(0))/dble(nxgrid))
      pdy=((yrange(1)-yrange(0))/dble(nygrid))      
!      hmpx=hmp*pdx
!      hmpy=hmp*pdy

! set grid data
      allocate(fval(0:nxgrid-1,0:nygrid-1,0:nsmapval))
   
      if(myrank.eq.0) then
        write(6,*) ' pixel size in x- and y-direction=',pdx,pdy
!        write(6,*) ' smoothing length in x- and y-direction=',hmpx,hmpy
      endif

! set how to read the data
      if(flagsp.eq.0.or.flagsp.eq.2) then
        flagr=1
      else
        flagr=0
      endif

      if(nstep.gt.1.and.myrank.eq.0) then
        open(49,file='./ini/file.dat',status='old')
      endif
      do istep=0,nstep-1
        if(nstep.gt.1.and.myrank.eq.0) then
          read(49,*) step
        endif    
        call MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! *** read the data ***
        call rdata(step,ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns,ai,tu,flagr)
        tu=tu*TMUGYR
        if(myrank.eq.0) then
           write(6,*) ' step, t(GYR)=',step,tu
           write(6,*) ' ngt,ndmt,ndm1t,nst=',ngt,ndmt,ndm1t,nst
        endif
        np=ng+ns+1

        if(flagrot.ne.0) then
! change rotation axis
          if(flagsp.eq.0.or.flagsp.eq.1) then
            do i=0,np-1
              y_p(i)=-y_p(i)
              vy_p(i)=-vy_p(i)
              vz_p(i)=-vz_p(i)
            enddo
          else
            do i=0,ndm-1
              y_dm(i)=-y_dm(i)
              vy_dm(i)=-vy_dm(i)
              vz_dm(i)=-vz_dm(i)
            enddo
          endif
        endif

        if(flag3d.lt.0) then
          if(flagsp.eq.0.or.flagsp.eq.1) then
            allocate(ty(0:np-1))
            allocate(tz(0:np-1))
! switch y and z
            do i=0,np-1
              ty(i)=y_p(i)
              tz(i)=z_p(i)
              y_p(i)=tz(i)
              z_p(i)=ty(i)
! for velocity
              ty(i)=vy_p(i)
              tz(i)=vz_p(i)
              vy_p(i)=tz(i)
              vz_p(i)=ty(i)
            enddo
            deallocate(ty)
            deallocate(tz)
          else
            allocate(ty(0:ndm-1))
            allocate(tz(0:ndm-1))
! switch y and z
            do i=0,ndm-1
              ty(i)=y_dm(i)
              tz(i)=z_dm(i)
              y_dm(i)=tz(i)
              z_dm(i)=ty(i)
! for velocity
              ty(i)=vy_dm(i)
              tz(i)=vz_dm(i)
              vy_dm(i)=tz(i)
              vz_dm(i)=ty(i)
            enddo
            deallocate(ty)
            deallocate(tz)
          endif
        endif             


        if(flagsp.eq.0) then
! *** gas
          nsp=0

          allocate(slist(0:ng-1))

          allocate(nh_p(0:np-1))
          allocate(temp_p(0:np-1))

! check h_p limit and set density and temparature, change units for position
          open(60,file='gastest.dat',status='unknown')
          do i=0,ng-1
            pn=list_ap(i)
            x_p(pn)=x_p(pn)*LUKPC
            y_p(pn)=y_p(pn)*LUKPC
            z_p(pn)=z_p(pn)*LUKPC
            h_p(pn)=h_p(pn)*LUKPC
            if(h_p(pn).lt.hmp) then
              h_p(pn)=hmp
            endif
! density and temperature
            p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn) 
            nh_p(pn)=((m_p(pn)-(mzZ_p(pn)+mzHe_p(pn))/MUSM) &
                       /m_p(pn))*rho_p(pn)*(DU/MP)
            temp_p(pn)=(p_p(pn)*myu_p(pn)/(rho_p(pn)*TPRHO*MYU))*1.0e4
            if(mod(i,10).eq.0) then
              write(60,'(2(1pE13.5))') nh_p(pn),temp_p(pn)
            endif
          enddo
          close(60)
 
          if(flag3d.le.0) then
            do i=0,ng-1
              pn=list_ap(i)
! make a list of particles around the region of the interest.
              if(flagid.eq.0.or.(id_p(pn).ge.ids.and.id_p(pn).le.ide)) then
!                if(nh_p(pn).gt.1.0.and.temp_p(pn).lt.1.0e4) then
                if(z_p(pn).gt.zsl-h_p(pn).and.z_p(pn).lt.zsl+h_p(pn)) then
                  if((x_p(pn).gt.xrange(0)-h_p(pn) &
                    .and.x_p(pn).lt.xrange(1)+h_p(pn)) &
                    .and.(y_p(pn).gt.yrange(0)-h_p(pn) &
                    .and.y_p(pn).lt.yrange(1)+h_p(pn))) then
                    slist(nsp)=pn
                    nsp=nsp+1
                  endif
                endif
!                endif
              endif
            enddo
          else
! 2D case, selected within the selected z region
            do i=0,ng-1
              pn=list_ap(i)
! make a list of particles around the region of the interest.
              if(flagid.eq.0.or.(id_p(pn).ge.ids.and.id_p(pn).le.ide)) then
!                if(nh_p(pn).gt.1.0.and.temp_p(pn).lt.1.0e4) then
                if(z_p(pn).gt.zrange(0).and.z_p(pn).lt.zrange(1)) then
                  if((x_p(pn).gt.xrange(0)-h_p(pn) &
                    .and.x_p(pn).lt.xrange(1)+h_p(pn)) &
                    .and.(y_p(pn).gt.yrange(0)-h_p(pn) &
                    .and.y_p(pn).lt.yrange(1)+h_p(pn))) then
                    z_p(pn)=0.0
                    slist(nsp)=pn
                    nsp=nsp+1
                  endif
                endif
!                endif
              endif
            enddo
          endif
          if(myrank.eq.0) then
            write(6,*) ' N seletcted gas particles=',nsp
          endif

          call allocate_smap_particle(nsp)

          open(60,file='test.dat',status='unknown')
          do i=0,nsp-1  
            pn=slist(i)
            xsp(i)=x_p(pn)
            ysp(i)=y_p(pn)
            zsp(i)=z_p(pn)
            vxsp(i)=vx_p(pn)*VUKMS
            vysp(i)=vy_p(pn)*VUKMS
            vzsp(i)=vz_p(pn)*VUKMS
! unit solar mass
            masssp(i)=m_p(pn)*MUSM
            psp(i)=p_p(pn)
            tempsp(i)=(p_p(pn)*myu_p(pn)/(rho_p(pn)*TPRHO*MYU))*1.0e4
            mhsp(i)=m_p(pn)*MUSM-(mzZ_p(pn)+mzHe_p(pn))
            metsp(i)=mzZ_p(pn)
            hsp(i)=h_p(pn)
            write(60,'(12(1pE13.5))') xsp(i),ysp(i),zsp(i) &
              ,vxsp(i),vysp(i),vzsp(i),masssp(i),psp(i),tempsp(i) &
              ,mhsp(i),metsp(i),hsp(i)
          enddo
          close(60)
          if(myrank.eq.0) then
            write(6,*) ' rank 0 N selected particles=',nsp
          endif

          deallocate(slist)

        endif

! initialisation
        do k=0,nsmapval-1
          do j=0,nygrid-1
            do i=0,nxgrid-1
              fval(i,j,k)=0.0d0 
            enddo
          enddo
        enddo

        call setsvalslice(nsp,nxgrid,nygrid,xrange(0),xrange(1) &
          ,yrange(0),yrange(1),nsmapval,flag3d,zsl,fval)

! covert the units
        do j=0,nygrid-1
          do i=0,nxgrid-1
            if(fval(i,j,0).gt.0.0d0) then
! ***  0:m, 1:vx, 2:vy, 3:vz, 4:p, 5:T, 6:mH, 7:Z
              do k=1,5
                fval(i,j,k)=fval(i,j,k)/fval(i,j,0)
              enddo
! nhp Msun/kpc^3
              k=6
              if(flag3d.gt.0) then
                fval(i,j,k)=fval(i,j,k)/(zrange(1)-zrange(2))
              endif
! mass weighted metallicity
              k=7
              fval(i,j,k)=fval(i,j,k)/fval(i,j,0)
            endif
          enddo
        enddo

! output
        if(myrank.eq.0) then
          if(flagsp.eq.0) then
            write(filename,'(a13,i6.6,a4)') 'output/slmapg',step,'.dat'
          else
            write(filename,'(a12,i6.6,a4)') 'output/slmap',step,'.dat'
          endif
          open(60,file=filename,status='unknown',form='unformatted')
          write(60) nxgrid,nygrid
          write(60) flag3d,flagrot
          write(60) xrange(0),xrange(1)
          write(60) yrange(0),yrange(1)
          write(60) zrange(0),zrange(1)
          write(60) nsmapval
          write(60) zsl
          write(60) tu
          write(60) hmp
          write(60) flagid
          write(60) ids,ide
! ascii output
          do k=0,nsmapval-1
            do j=0,nygrid-1
              write(60) (fval(i,j,k),i=0,nxgrid-1)
            enddo
          enddo
          close(60)
          write(filename,'(a12,i6.6,a4)') 'output/slmap',step,'.asc'
          open(61,file=filename,status='unknown')
          do j=0,nygrid-1
            ypm=yrange(0)+(dble(j)+0.5d0)*pdy
            do i=0,nxgrid-1
              xpm=xrange(0)+(dble(i)+0.5d0)*pdx
              write(61,'(10(1pE13.5))') xpm,ypm,fval(i,j,0),fval(i,j,1) &
                ,fval(i,j,2),fval(i,j,3),fval(i,j,4),fval(i,j,5) &
                ,fval(i,j,6),fval(i,j,7)
            enddo
          enddo
        endif 

      enddo

      deallocate(fval)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program
