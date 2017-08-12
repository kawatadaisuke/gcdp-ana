! /***************************************
!   smrvrot ver.1
!  12 Aug. 2017  written by D.Kawata
! ****************************************/

program smrvrot
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
   
      implicit none
      include 'mpif.h'

      integer id,i,j,k,pn,flag,ns,ng,np,istep
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer nstep,step,is,ic
      integer nval
      integer flagr,flagc
      integer ierr
! particle data
      double precision,allocatable :: rxyp(:),vrotp(:),agep(:)
      double precision,allocatable :: xsp(:),ysp(:)
! grid info
      integer nxgrid,nygrid
      double precision xrange(0:1),yrange(0:1)
      double precision hmp,hmpx,hmpy,pdx,pdy,xpm,ypm,xlsm,xhsm,ylsm,yhsm
      double precision,allocatable :: fval(:,:)
! gas properties
      double precision ltempp,lnhp
! gas properties selection
      integer flagselg
      double precision ltemprange(0:1),lnhrange(0:1)
! for stars
      integer flagsout
      double precision agerange(0:1)
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
        write(6,*) ' smapR_Vrot ver. 1 12/08/17'
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
        read(50,*) flagselg
        if(flagselg.ne.0) then
          read(50,*) ltemprange(0),ltemprange(1)
          read(50,*) lnhrange(0),lnhrange(1)
        else
          ltemprange(0)=-INF
          ltemprange(1)=INF
          lnhrange(0)=-INF
          lnhrange(1)=INF
        endif
        read(50,*) flagsout
        if(flagsout.ne.0) then
          read(50,*) agerange(0),agerange(1)
        endif
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
        write(6,*) ' smoothing length if not gas (pixel)=',hmp
        if(flagselg.ne.0) then
          write(6,*) ' Log T range=',ltemprange(0),ltemprange(1)
          write(6,*) ' Log nh range=',lnhrange(0),lnhrange(1)
        endif
        if(flagsout.ne.0) then
          write(6,*) ' output stellar R and Vrot for stars age range =' &
            ,agerange(0),agerange(1)
        endif
! send the data to the other node
      endif

      nval=6
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        tivr(2)=nxgrid
        tivr(3)=nygrid
        tivr(4)=flagselg
        tivr(5)=flagsout
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      nxgrid=tivr(2)
      nygrid=tivr(3)
      flagselg=tivr(4)
      flagsout=tivr(5)
      deallocate(tivr)

      nval=11
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=xrange(0)
        tdvr(1)=xrange(1)
        tdvr(2)=yrange(0)
        tdvr(3)=yrange(1)
        tdvr(4)=hmp
        tdvr(5)=ltemprange(0)
        tdvr(6)=ltemprange(1)
        tdvr(7)=lnhrange(0)
        tdvr(8)=lnhrange(1)
        tdvr(9)=agerange(0)
        tdvr(10)=agerange(1)
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      xrange(0)=tdvr(0)
      xrange(1)=tdvr(1)
      yrange(0)=tdvr(2)
      yrange(1)=tdvr(3)
      hmp=tdvr(4)
      ltemprange(0)=tdvr(5)
      ltemprange(1)=tdvr(6)
      lnhrange(0)=tdvr(7)
      lnhrange(1)=tdvr(8)
      agerange(0)=tdvr(9)
      agerange(1)=tdvr(10)

      deallocate(tdvr)

! smoothing length for x and y
      pdx=((xrange(1)-xrange(0))/dble(nxgrid))
      pdy=((yrange(1)-yrange(0))/dble(nygrid))      
      hmpx=hmp*pdx
      hmpy=hmp*pdy
! set grid data
      allocate(fval(0:nxgrid-1,0:nygrid-1))
   
      if(myrank.eq.0) then
        write(6,*) ' pixel size in x- and y-direction=',pdx,pdy
        write(6,*) ' smoothing length in x- and y-direction=',hmpx,hmpy
      endif

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

        if(flagc.eq.0) then
! *** gas
          np=0

          allocate(slist(0:ng-1))
          allocate(vrotp(0:ng+ns-1))
          allocate(rxyp(0:ng+ns-1))

        open(60,file='test.dat',status='unknown')
! make a list of particles around the region of the interest.
          do i=0,ng-1
            pn=list_ap(i)
            rxyp(pn)=dsqrt(x_p(pn)**2+y_p(pn)**2)*LUKPC
            if(rxyp(pn).lt.xrange(1)+hmpx) then
              if(rxyp(pn).gt.0.0d0) then
                vrotp(pn)=((vx_p(pn)*y_p(pn)-vy_p(pn)*x_p(pn)) &
                  /(rxyp(pn)/LUKPC))*VUKMS
              else
                vrotp(pn)=0.0d0
              endif
              if(vrotp(pn).gt.yrange(0)-hmpy.and. &
                 vrotp(pn).lt.yrange(1)+hmpy) then
                if(flagselg.ne.0) then
                  lnhp=dlog10(((m_p(pn)-((mzZ_p(pn)+mzHe_p(pn))/MUSM)) &
                       /m_p(pn))*rho_p(pn)*(DU/MP))
                  p_p(pn)=(GAM-1.0d0)*rho_p(pn)*u_p(pn) 
                  ltempp=dlog10((p_p(pn)*myu_p(pn)/(rho_p(pn)*TPRHO*MYU))*1.0e4)
                  if((lnhp.gt.lnhrange(0).and.lnhp.lt.lnhrange(1)).and. &
                     (ltempp.gt.ltemprange(0).and.ltempp.lt.ltemprange(1))) then
                    slist(np)=pn
                    np=np+1

                    write(60,'(7(1pE13.5))') x_p(pn),y_p(pn),z_p(pn) &
                      ,rxyp(pn),vrotp(pn),lnhp,ltempp

                  endif
                else
                  slist(np)=pn
                  np=np+1
                endif
              endif
            endif
          enddo
          close(60)
          if(myrank.eq.0) then
            write(6,*) ' rank 0 N selected particles=',np
          endif
        endif

        allocate(xsp(0:np-1))
        allocate(ysp(0:np-1))

! set particle for smoothed map
        do i=0,np-1
          pn=slist(i)
          xsp(i)=(rxyp(pn)-xrange(0))/pdx
          ysp(i)=(vrotp(pn)-yrange(0))/pdy
        enddo

        deallocate(slist)

! initialisation
        do j=0,nygrid-1
          do i=0,nxgrid-1
            fval(i,j)=0.0d0 
          enddo
        enddo

! set smoothing
        xlsm=0.0d0
        xhsm=dble(nxgrid)
        ylsm=0.0d0
        yhsm=dble(nygrid)

        call setsval(np,xsp,ysp,hmp,nxgrid,nygrid,xlsm,xhsm,ylsm,yhsm,fval)

! output
        if(myrank.eq.0) then
          open(60,file='output/smapR_Vrot.dat',status='unknown' &
              ,form='unformatted')
          open(61,file='output/smapR_Vrot.asc',status='unknown')
          write(60) nxgrid,nygrid
          write(60) xrange(0),xrange(1)
          write(60) yrange(0),yrange(1)
          do j=0,nygrid-1
            ypm=yrange(0)+(dble(j)+0.5d0)*pdy
            do i=0,nxgrid-1
              xpm=xrange(0)+(dble(i)+0.5d0)*pdx
              write(60) xpm,ypm,fval(i,j)
              write(61,'(3(1pE13.5))') xpm,ypm,fval(i,j)
            enddo
          enddo
          close(60)
        endif 

        deallocate(rxyp)
        deallocate(vrotp)
        deallocate(xsp)
        deallocate(ysp)

! star R vrot data output
        if(flagsout.ne.0.and.myrank.eq.0) then
          open(60,file='output/starR_Vrot.asc',status='unknown')

          allocate(agep(0:ng+ns-1))
          allocate(rxyp(0:ng+ns-1))
          allocate(vrotp(0:ng+ns-1))
          allocate(slist(0:ns-1))

! make a list of particles around the region of the interest.
          np=0
          do i=ng,ng+ns-1
            pn=list_ap(i)
            agep(pn)=(tu-ts_p(pn))*TMUGYR
            if(agep(pn).ge.agerange(0).and.agep(pn).le.agerange(1)) then
              slist(np)=pn
              np=np+1
            endif
          enddo

          write(60,'(a5,I10)') '# np=',np
!                            12345678901234567890
          write(60,'(a20)') '# x y z Rxy Vrot age'
          do i=0,np-1
            pn=slist(i)
            rxyp(pn)=dsqrt(x_p(pn)**2+y_p(pn)**2)*LUKPC
            if(rxyp(pn).gt.0) then
              vrotp(pn)=((vx_p(pn)*y_p(pn)-vy_p(pn)*x_p(pn)) &
                /(rxyp(pn)/LUKPC))*VUKMS
            else
              vrotp(pn)=0.0d0
            endif
            write(60,'(6(1pE13.5))') x_p(pn),y_p(pn),z_p(pn) &
                ,rxyp(pn),vrotp(pn),agep(pn)
          enddo
          close(60)
        endif
      enddo

      deallocate(fval)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program
