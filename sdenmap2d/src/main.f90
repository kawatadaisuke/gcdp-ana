! /***************************************
!   sdenmap2d ver.9.1
!  1 Aug. 2017  written by D.Kawata
! ****************************************/

program sdenmap2d
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
      use particle
      use mesh
   
      implicit none
      include 'mpif.h'

      integer id,i,j,k,pn,flag,ns,ng,np,istep
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer nstep,step,is,ic
      integer nval
      integer flagc,flagcom,flagr,flagidr
      integer idr(0:1)
      integer ierr
      integer flagrot
      double precision degrot
! inclination
      integer flaginc
      double precision iang
! for box grid info
      integer nxgrid,nygrid
      double precision xrange(0:1),yrange(0:1)
      double precision hmp,dpp,angp,hmin,rpi
! galb selection
      integer flaggalbr
      double precision galbrange(0:1),galbp,rsun,degsun,rxyp
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! for adjusting COM
      integer flagcmx,nsp,nc,ncr,nct
      double precision cx,cy,cz,mt,rmax0,cxt,cyt,czt,rmax
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
        write(6,*) ' sdenmap2d ver. p9 28/07/16'
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
        read(50,*) nxgrid,nygrid
        read(50,*) xrange(0),xrange(1)
        read(50,*) yrange(0),yrange(1)
        read(50,*) flagc
        read(50,*) hmp
        read(50,*) flagcmx
        read(50,*) rmax0   
        read(50,*) flaginc
        read(50,*) iang
        read(50,*) flagidr
        read(50,*) idr(0),idr(1)
        if(flagidr.eq.0) then
          idr(0)=0
          idr(1)=1000000000
        endif
        read(50,*) flaggalbr
        read(50,*) galbrange(0),galbrange(1)
        read(50,*) rsun,degsun
        read(50,*) flagrot
        read(50,*) degrot
        close(50)
        if(nstep.gt.1) then
          write(6,*) ' nstep=',nstep
        else
          write(6,*) ' step=',step 
        endif
        write(6,*) ' grid nx, ny=',nxgrid,nygrid
        write(6,*) ' x range (kpc)=',xrange(0),xrange(1)
        write(6,*) ' y range (kpc)=',yrange(0),yrange(1)
        if(flagc.eq.0) then
          write(6,*) ' for gas'
        else if(flagc.eq.1) then
          write(6,*) ' for dm'
        else
          write(6,*) ' for star'
        endif
        write(6,*) ' smoothing length if not gas (kpc)=',hmp
        if(flagcmx.ne.0) then
          if(flagcmx.eq.1) then
            write(6,*) ' adjust centre of the mass with baryon particles.'
          else
            write(6,*) ' adjust centre of the mass with DM particles.'
          endif
          write(6,*) ' initial rmax=',rmax0
        endif
        if(flaginc.ne.0) then
          write(6,*) ' rotate the galaxy along x-axis by iang=',iang
        endif
        if(flagidr.ne.0) then
          write(6,*) ' use the particles only in the ID range between' &
           ,idr(0),' and ',idr(1),' (kpc)'
        endif
        if(flaggalbr.ne.0) then
          write(6,*) ' only particle within the brange between=' &
           ,galbrange(0),galbrange(1)
          write(6,*) ' rotate degsun=',degsun
          write(6,*) ' put the observer at (-',rsun,',0) (kpc)'
        endif
        if(flagrot.ne.0) then
          if(flagrot.lt.0) then
            write(6,*) ' rotation axis set to be -z'
          endif
          write(6,*) ' rotate ',degrot,' (deg)'
       endif
       if(flaggalbr.ne.0.and.flagrot.ne.0) then
         write(6,*) ' Error: not work with flaggalbr.ne.0 and flagrot.ne.0.'
         stop
       endif
! send the data to the other node
      endif

      nval=12
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        tivr(2)=nxgrid
        tivr(3)=nygrid
        tivr(4)=flagc
        tivr(5)=flagidr
        tivr(6)=idr(0)
        tivr(7)=idr(1)
        tivr(8)=flagcmx
        tivr(9)=flaginc
        tivr(10)=flaggalbr
        tivr(11)=flagrot
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      nxgrid=tivr(2)
      nygrid=tivr(3)
      flagc=tivr(4)
      flagidr=tivr(5)
      idr(0)=tivr(6)
      idr(1)=tivr(7)
      flagcmx=tivr(8)
      flaginc=tivr(9)
      flaggalbr=tivr(10)
      flagrot=tivr(11)
      deallocate(tivr)

      nval=12
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=xrange(0)
        tdvr(1)=xrange(1)
        tdvr(2)=yrange(0)
        tdvr(3)=yrange(1)
        tdvr(4)=hmp
        tdvr(5)=rmax0
        tdvr(6)=iang
        tdvr(7)=galbrange(0)
        tdvr(8)=galbrange(1)
        tdvr(9)=rsun
        tdvr(10)=degsun
        tdvr(11)=degrot
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      xrange(0)=tdvr(0)
      xrange(1)=tdvr(1)
      yrange(0)=tdvr(2)
      yrange(1)=tdvr(3)
      hmp=tdvr(4)
      rmax0=tdvr(5)
      iang=tdvr(6)
      galbrange(0)=tdvr(7)
      galbrange(1)=tdvr(8)
      rsun=tdvr(9)
      degsun=tdvr(10)
      degrot=tdvr(11)

      deallocate(tdvr)

! set how to read the data
      if(flagc.eq.0.or.flagc.eq.2) then
        flagr=1
      else
        flagr=0
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
        if(myrank.eq.0) then
           write(6,*) ' step, t(GYR)=',step,tu*TMUGYR
           write(6,*) ' ngt,ndmt,ndm1t,nst=',ngt,ndmt,ndm1t,nst
        endif
        if(flagidr.eq.0.and.idr(1).lt.ngt+nst) then
          idr(1)=2*(ngt+nst)
        endif
        if(flagidr.eq.0.and.idr(1).lt.ndmt) then
          idr(1)=2*ndmt
        endif

        np=ng+ns
! adjust the centre to the centre of the mass
        if(flagcmx.ne.0) then
! minimum number of particles to determine COM
          ncr=50     
          nsp=0
          nval=4
          allocate(tdvr(0:nval-1))
          allocate(tdvs(0:nval-1))
          if(flagcmx.eq.1) then
            allocate(tx(0:np-1))
            allocate(ty(0:np-1))
            allocate(tz(0:np-1))
            allocate(tm(0:np-1))
            allocate(rp(0:np-1))
            do i=0,np-1
              pn=list_ap(i)
              tx(nsp)=x_p(pn)
              ty(nsp)=y_p(pn)
              tz(nsp)=z_p(pn)
              tm(nsp)=m_p(pn)
              nsp=nsp+1
            enddo
          else
            allocate(tx(0:ndm-1))
            allocate(ty(0:ndm-1))
            allocate(tz(0:ndm-1))
            allocate(tm(0:ndm-1))
            allocate(rp(0:ndm-1))
            do i=0,ndm-1
              tx(nsp)=x_dm(i)
              ty(nsp)=y_dm(i)
              tz(nsp)=z_dm(i)
              tm(nsp)=m_dm(i)   
              nsp=nsp+1
            enddo
          endif
          if(myrank.eq.0) then
            write(6,*) ' initial nsp for adjusting COM in myrank 0=',nsp
          endif
          cx=0.0d0
          cy=0.0d0
          cz=0.0d0
          rmax=rmax0
          do j=1,90
            nc=0 
            cxt = 0.0d0
            cyt = 0.0d0
            czt = 0.0d0
            mt = 0.0d0
            do i =0,nsp-1
              rp(i)=dsqrt(tx(i)*tx(i)+ty(i)*ty(i)+tz(i)*tz(i))
              if(rp(i).lt.rmax) then
                cxt=cxt+tm(i)*tx(i)
                cyt=cyt+tm(i)*ty(i)
                czt=czt+tm(i)*tz(i)
                mt=mt+tm(i)
                nc=nc+1
              endif
            enddo
            nct=0
            call MPI_ALLREDUCE(nc,nct,1,MPI_INTEGER,MPI_SUM &
             ,MPI_COMM_WORLD,ierr)
            nc=nct
            if(nct.le.ncr) then
              goto 90
            endif
            tdvs(0)=cxt
            tdvs(1)=cyt
            tdvs(2)=czt
            tdvs(3)=mt
            do i=0,nval-1  
              tdvr(i)=0.0d0
            enddo
            call MPI_ALLREDUCE(tdvs,tdvr,nval,MPI_DOUBLE_PRECISION &
              ,MPI_SUM,MPI_COMM_WORLD,ierr)
            cxt=tdvr(0)
            cyt=tdvr(1)
            czt=tdvr(2)
            mt=tdvr(3)

            cxt=cxt/mt
            cyt=cyt/mt
            czt=czt/mt
            do i =0,nsp-1
              tx(i)=tx(i)-cxt
              ty(i)=ty(i)-cyt 
              tz(i)=tz(i)-czt
            enddo
            cx=cx+cxt
            cy=cy+cyt
            cz=cz+czt
!            if(myrank.eq.0) then
!              write(6,*) j,' rmax,nc,cx=',rmax,nc,cx,cy,cz
!            endif
            if(mod(j,3).eq.0) then
              rmax = rmax*0.9
            endif
          enddo
 90       if(myrank.eq.0) then
            write(6,*) ' rmax,nc,cx=',rmax,nc,cx,cy,cz
          endif
          do i=0,np-1
            x_p(i)=x_p(i)-cx
            y_p(i)=y_p(i)-cy
            z_p(i)=z_p(i)-cz
          enddo
!!          open(60,file='tdmp.dat', status='unknown')
          do i=0,ndm-1
            x_dm(i)=x_dm(i)-cx
            y_dm(i)=y_dm(i)-cy
            z_dm(i)=z_dm(i)-cz
!            if(mod(i,10).eq.0) then
!              write(60,'(3(1pE13.5))') x_dm(i),y_dm(i),z_dm(i)
!            endif
          enddo
!          close(60)

          deallocate(tdvs)
          deallocate(tdvr)
          deallocate(tx)
          deallocate(ty)
          deallocate(tz)
          deallocate(tm)
          deallocate(rp)
        endif

! rotate the system
        if(flaggalbr.ne.0) then
          if(myrank.eq.0) then
            write(6,*) ' rotate the system by ',degsun,' (deg)'
          endif
! deg -> radian
          degsun=M_PI*degsun/180.0d0
! for baryon
          allocate(tx(0:np-1))
          allocate(ty(0:np-1))

          do i=0,np-1
            tx(i)=x_p(i)*dcos(degsun)-y_p(i)*dsin(degsun)
            ty(i)=x_p(i)*dsin(degsun)+y_p(i)*dcos(degsun)
            x_p(i)=tx(i)
            y_p(i)=ty(i)
          enddo
          close(60)

! for DM
          deallocate(tx)
          deallocate(ty)
          allocate(tx(0:ndm-1))
          allocate(ty(0:ndm-1))

!          open(60,file='rotdm.dat',status='unknown')
          do i=0,ndm-1
            tx(i)=x_dm(i)*dcos(degsun)-y_dm(i)*dsin(degsun)
            ty(i)=x_dm(i)*dsin(degsun)+y_dm(i)*dcos(degsun)
            x_dm(i)=tx(i)
            y_dm(i)=ty(i)

!            write(60,'(2(1pE13.5))') x_dm(i),y_dm(i)

          enddo
!          close(60)

          deallocate(tx)
          deallocate(ty)

        endif

! position kpc and velocity
        do i=0,np-1
          x_p(i)=x_p(i)*LUKPC
          y_p(i)=y_p(i)*LUKPC
          z_p(i)=z_p(i)*LUKPC
          vx_p(i)=vx_p(i)*VUKMS
          vy_p(i)=vy_p(i)*VUKMS
          vz_p(i)=vz_p(i)*VUKMS
! smoothing length
          h_p(i)=h_p(i)*LUKPC
        enddo
        do i=0,ndm-1
          x_dm(i)=x_dm(i)*LUKPC
          y_dm(i)=y_dm(i)*LUKPC
          z_dm(i)=z_dm(i)*LUKPC
          vx_dm(i)=vx_dm(i)*VUKMS
          vy_dm(i)=vy_dm(i)*VUKMS
          vz_dm(i)=vz_dm(i)*VUKMS
        enddo

        if(flagc.eq.0) then
! *** gas
          np=0

          allocate(slist(0:ng-1))

!        open(60,file='test.dat',status='unknown')
! make a list of particles around the region of the interest.
          do i=0,ng-1
            pn=list_ap(i)
            if(x_p(pn).gt.1.5d0*xrange(0)-hmp.and.x_p(pn).lt.1.5d0*xrange(1)+hmp) then
            if(y_p(pn).gt.1.5d0*yrange(0)-hmp.and.y_p(pn).lt.1.5d0*yrange(1)+hmp) then
              if(id_p(pn).ge.idr(0).and.id_p(pn).le.idr(1)) then
              if(flagfd_p(pn).eq.0) then
              if(flaggalbr.eq.0) then
                slist(np)=pn
                np=np+1
              else
                rxyp=dsqrt((x_p(pn)+rsun)**2+y_p(pn)**2)
                if(rxyp.gt.0.0d0) then
                  galbp=datan(z_p(pn)/rxyp)*180.0d0/M_PI
                else
                  galbp=0.0d0
                endif
                if(galbp.gt.galbrange(0).and.galbp.lt.galbrange(1)) then
                  slist(np)=pn
                  np=np+1
                endif
              endif
              endif
              endif
            endif
            endif
          enddo

          call allocate_particle(np)

          do i=0,np-1
            pn=slist(i)
            xp(i)=x_p(pn)
            yp(i)=y_p(pn)
            zp(i)=z_p(pn)
            vxp(i)=vx_p(pn)
            vyp(i)=vy_p(pn)
            vzp(i)=vz_p(pn)
            massp(i)=m_p(pn)
            metp(i)=mzZ_p(pn)/(m_p(pn)*MUSM)
            hp(i)=h_p(pn)
            idp(i)=id_p(pn)
          enddo

          deallocate(slist)

        else if(flagc.eq.1) then
! *** DM
          allocate(slist(0:ndm1-1))

!          open(60,file='dmgalb.dat',status='unknown')
! make a list of particles around the region of the interest.
          np=0
          do i=0,ndm1-1
            pn=list_adm(i)
            if(x_dm(pn).gt.1.5d0*xrange(0)-hmp.and.x_dm(pn).lt.1.5d0*xrange(1)+hmp) then
            if(y_dm(pn).gt.1.5d0*yrange(0)-hmp.and.y_dm(pn).lt.1.5d0*yrange(1)+hmp) then
              if(id_dm(pn).ge.idr(0).and.id_dm(pn).le.idr(1)) then
              if(flaggalbr.eq.0) then
                slist(np)=pn
                np=np+1
              else
                rxyp=dsqrt((x_dm(pn)+rsun)**2+y_dm(pn)**2)
                if(rxyp.gt.0.0d0) then
                  galbp=datan(z_dm(pn)/rxyp)*180.0d0/M_PI
                else
                  galbp=0.0d0
                endif
                if(galbp.gt.galbrange(0).and.galbp.lt.galbrange(1)) then
                  slist(np)=pn
                  np=np+1
                endif

!                write(60,'(4(1pE13.5))') x_dm(pn),y_dm(pn),z_dm(pn),galbp
 
              endif
              endif
            endif
            endif
          enddo

!          close(60)

          call allocate_particle(np)

          do i=0,np-1
            pn=slist(i)
            xp(i)=x_dm(pn)
            yp(i)=y_dm(pn)
            zp(i)=z_dm(pn)
            vxp(i)=vx_dm(pn)
            vyp(i)=vy_dm(pn)
            vzp(i)=vz_dm(pn)
            massp(i)=m_dm(pn)
            hp(i)=h_dm(pn)
            metp(i)=0.0d0
            idp(i)=id_dm(pn)
          enddo

          deallocate(slist)

        else
! *** star
          allocate(slist(0:ns-1))

! make a list of particles around the region of the interest.
          np=0
          do i=0,ns-1
            pn=list_ap(i+ng)
            if(x_p(pn).gt.1.5d0*xrange(0)-hmp.and.x_p(pn).lt.1.5d0*xrange(1)+hmp) then
            if(y_p(pn).gt.1.5d0*yrange(0)-hmp.and.y_p(pn).lt.1.5d0*yrange(1)+hmp) then
              if(id_p(pn).ge.idr(0).and.id_p(pn).le.idr(1)) then
              if(flaggalbr.eq.0) then
                slist(np)=pn
                np=np+1
              else
                rxyp=dsqrt((x_p(pn)+rsun)**2+y_p(pn)**2)
                if(rxyp.gt.0.0d0) then
                  galbp=datan(z_p(pn)/rxyp)*180.0d0/M_PI
                else
                  galbp=0.0d0
                endif
                if(galbp.gt.galbrange(0).and.galbp.lt.galbrange(1)) then
                  slist(np)=pn
                  np=np+1
                endif
              endif
              endif
            endif
            endif
          enddo

          call allocate_particle(np)

          do i=0,np-1
            pn=slist(i)
            xp(i)=x_p(pn)
            yp(i)=y_p(pn)
            zp(i)=z_p(pn)
            vxp(i)=vx_p(pn)
            vyp(i)=vy_p(pn)
            vzp(i)=vz_p(pn)
            massp(i)=m_p(pn)
            hp(i)=((XHSOL*(DU/MP)/NSTH)**(1.0d0/3.0d0)*ETAH)*(m_p(pn)**THIRD)
            metp(i)=mzZ_p(pn)/(m_p(pn)*MUSM)
            idp(i)=id_p(pn)
          enddo

          deallocate(slist)
        endif
! getting the total np
        call MPI_ALLREDUCE(np,npt,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

        if(myrank.eq.0) then   
          write(6,*) 'Selected npt, np(myrank)=',npt,np
        endif

        if(flagrot.ne.0) then
! rotate the system
          if(flagrot.lt.0) then
            do i=0,np-1
              yp(i)=-yp(i)
              zp(i)=-zp(i)
              vyp(i)=-vyp(i)
              vzp(i)=-vzp(i)
            enddo
          endif

! convert deg -> radian
          degrot=degrot*M_PI/180.0d0
! rotate the galaxy along x-axis
          allocate(tx(0:np-1))
          allocate(ty(0:np-1))

          do i=0,np-1
! position
            tx(i)=xp(i)*dcos(-degrot)-yp(i)*dsin(-degrot)
            ty(i)=xp(i)*dsin(-degrot)+yp(i)*dcos(-degrot)
            xp(i)=tx(i)
            yp(i)=ty(i)
! velocity
            tx(i)=vxp(i)*dcos(-degrot)-vyp(i)*dsin(-degrot)
            ty(i)=vxp(i)*dsin(-degrot)+vyp(i)*dcos(-degrot)
            vxp(i)=tx(i)
            vyp(i)=ty(i)
          enddo

          deallocate(tx)
          deallocate(ty)

        endif
    
! calculate vrot and vrad
        do i=0,np-1
          rpi=dsqrt(xp(i)**2+yp(i)**2)
          if(rpi.gt.0.0d0) then
            vradp(i)=(vxp(i)*xp(i)+vyp(i)*yp(i))/rpi
            vrotp(i)=(vxp(i)*yp(i)-vyp(i)*xp(i))/rpi
          else
            vradp(i)=0.0d0
            vrotp(i)=0.0d0
          endif
        enddo

        if(flaginc.ne.0) then

! convert deg -> radian
          iang=iang*M_PI/180.0d0
! rotate the galaxy along x-axis
          allocate(ty(0:np-1))
          allocate(tz(0:np-1))

          do i=0,np-1
! position
            ty(i)=yp(i)*dcos(iang)-zp(i)*dsin(iang)
            tz(i)=yp(i)*dsin(iang)+zp(i)*dcos(iang)
            yp(i)=ty(i)
            zp(i)=tz(i)
! velocity
            ty(i)=vyp(i)*dcos(iang)-vzp(i)*dsin(iang)
            tz(i)=vyp(i)*dsin(iang)+vzp(i)*dcos(iang)
            vyp(i)=ty(i)
            vzp(i)=tz(i)
          enddo

          deallocate(ty)
          deallocate(tz)

        endif


! no hp finding
!        if(flagc.ne.0) then
! use gcd+ gtree
!          call gtreebuild(np)
!          if(myrank.eq.0) then
!            write(6,*) ' start setvp'
!          endif
! use set_value to set hp(), ETAH to control the size of the smoothing length
!          call set_value(np)
!        endif

! set hmin 2x grid size
        hmin=2.0d0*(xrange(1)-xrange(0))/dble(nxgrid)
        if(hmin.lt.2.0d0*(yrange(1)-yrange(0))/dble(nygrid)) then
          hmin=2.0d0*(yrange(1)-yrange(0))/dble(nygrid)
        endif
        if(myrank.eq.0.and.istep.eq.1) then
           write(6,*) ' hmin=',hmin
        endif

! set hmp and hmax
        if(flagc.ne.0) then
          do i=0,np-1
            hp(i)=hmp
          enddo
        endif
        do i=0,np-1
          if(hp(i).lt.hmin) then
            hp(i)=hmin     
          endif
        enddo

! for test, not parallerize
!        write(fileo,'(a4,i3.3)') 'comp',myrank
!        open(60,file=fileo,status='unknown')
!        do i=0,np-1
!          write(60,'(8(1pE13.5),I10)') xp(i),yp(i),zp(i),vxp(i),vyp(i),vzp(i) &
!           ,hp(i),metp(i),idp(i)
!        enddo
!        close(60)

        call setsval(np,nxgrid,nygrid,xrange,yrange)

        call output(step,nxgrid,nygrid,xrange,yrange,tu,flagc,flaginc,iang)

! deallocate mesh data
        deallocate(x_m)
        deallocate(y_m)
        deallocate(denxy_m)
        deallocate(denxz_m)
        deallocate(vrot_m)
        deallocate(vrad_m)
        deallocate(vz_m)
        deallocate(met_m)

      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program sdenmap2d
