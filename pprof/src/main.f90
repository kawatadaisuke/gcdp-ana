! /***************************************
!   prof pver6
!  17 Apr. 2017  written by D.Kawata
! ****************************************/

program prof
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
      use particle
      use mesh
   
      implicit none
      include 'mpif.h'

      integer i,j,k,pn,idc,istep,idr
      integer flag,ns,ng,np
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer nstep,step
      integer nval,nvpr
      integer flagc,flagcom,flagr
      integer ierr
! id range for each component
      integer nidc
      integer nprs,nprst
      integer,allocatable :: idspc(:),ids(:),ide(:)
      double precision,allocatable :: agel(:),ageh(:)
! for profile
      integer ndr2d
      double precision rmin2d,rmax2d,dr,ri,ro,dvr
      double precision,allocatable :: dmr(:),vrm(:),vr2m(:) &
       ,vtm(:),vt2m(:),vzm(:),vz2m(:),amtr(:),z2m(:),mtr(:)
      double precision vrp,vtp,vtzp,rp
      double precision vrsig,vtsig,vzsig,rdmr,zhsig
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! ***  for filename ***
      character filename*60,fileo*60
! for work
      integer,allocatable :: slist(:)
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:),tdvs(:)

! ***** MPI Initialization *****
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      if(myrank.eq.0) then
        write(6,*) ' prof ver. p6 17/04/17'
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif
      if(nprocs.gt.NCPU) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: increase NCPU=',NCPU
        endif
        call MPI_FINALIZE()
        stop
      endif

! *****   Open Input File   ******
      if(myrank.eq.0) then
        open(50,file='./ini/input.dat',status='old')
        read(50,*) nstep,step
        read(50,*) ndr2d
        read(50,*) rmin2d,rmax2d
        read(50,*) nidc
        if(nstep.gt.1) then
          write(6,*) ' nstep=',nstep
        else
          write(6,*) ' step=',step 
        endif
        write(6,*) ' ndr,rmin,rmax=',ndr2d,rmin2d,rmax2d
        write(6,*) ' nidc=',nidc
      endif
      nval=4
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        tivr(2)=ndr2d
        tivr(3)=nidc
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      ndr2d=tivr(2)
      nidc=tivr(3)
      deallocate(tivr)
! allocate idrange
      allocate(idspc(0:nidc-1))
      allocate(ids(0:nidc-1))
      allocate(ide(0:nidc-1))
      allocate(agel(0:nidc-1))
      allocate(ageh(0:nidc-1))
      if(myrank.eq.0) then
        do idc=0,nidc-1
          read(50,*) idspc(idc)
          read(50,*) ids(idc),ide(idc)
          write(6,*) idc,'  idspc=',idspc(idc)
          write(6,*) ' id range=',ids(idc),ide(idc)
          if(idspc(idc).eq.2) then
            read(50,*) agel(idc),ageh(idc)
            write(6,*) ' age range for stars (Gyr) =',agel(idc),ageh(idc)
          endif
        enddo
        close(50)
      endif

! sending the id range information to the other procs
      nval=3*nidc
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        do idc=0,nidc-1
          tivr(idc*3)=idspc(idc)
          tivr(idc*3+1)=ids(idc)
          tivr(idc*3+2)=ide(idc)
        enddo
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      do idc=0,nidc-1
        idspc(idc)=tivr(idc*3)
        ids(idc)=tivr(idc*3+1)
        ide(idc)=tivr(idc*3+2)
      enddo
      deallocate(tivr)

      nval=2
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=rmin2d
        tdvr(1)=rmax2d
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rmin2d=tdvr(0)
      rmax2d=tdvr(1)
      deallocate(tdvr)

! for stellar age range
      nval=2*nidc
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        do idc=0,nidc-1
          tdvr(idc*2)=agel(idc)
          tdvr(idc*2+1)=ageh(idc)
        enddo
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      do idc=0,nidc-1
        agel(idc)=tdvr(idc*2)
        ageh(idc)=tdvr(idc*2+1)
      enddo
      deallocate(tdvr)

! calculating the radial range and bin size
      dr=(rmax2d-rmin2d)/dble(ndr2d)
! allocate radial profile bins
      allocate(dmr(0:ndr2d-1))
      allocate(vrm(0:ndr2d-1))
      allocate(vr2m(0:ndr2d-1))
      allocate(vtm(0:ndr2d-1))
      allocate(vt2m(0:ndr2d-1))
      allocate(vzm(0:ndr2d-1))
      allocate(vz2m(0:ndr2d-1))
      allocate(z2m(0:ndr2d-1))
      allocate(amtr(0:ndr2d-1))
      allocate(mtr(0:ndr2d-1))

! ste how to read the data
      if(flagc.eq.0) then
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

! analysing the profile for each component
        do idc=0,nidc-1
! initialise
          nprs=0
          do i=0,ndr2d-1
            dmr(i)=0.0d0
            vrm(i)=0.0d0
            vr2m(i)=0.0d0
            vtm(i)=0.0d0
            vt2m(i)=0.0d0
            vzm(i)=0.0d0
            vz2m(i)=0.0d0
            z2m(i)=0.0d0
            amtr(i)=0.0d0
            mtr(i)=0.0d0
          enddo
          if(idspc(idc).eq.0) then
! for gas
            do i=0,ng-1
              pn=list_ap(i)
              if(id_p(pn).ge.ids(idc).and.id_p(pn).le.ide(idc)) then
                rp=dsqrt(x_p(pn)**2+y_p(pn)**2)
                idr=int((rp-rmin2d)/dr)
                if(idr.ge.0.and.idr.lt.ndr2d) then
                  nprs=nprs+1
                  dmr(idr)=dmr(idr)+m_p(pn)
! calculate velocity components
                  if(rp.gt.0.0d0) then
                    vrp=(vx_p(pn)*x_p(pn)+vy_p(pn)*y_p(pn))/rp
! clock-wise
                    vtp=(vx_p(pn)*y_p(pn)-vy_p(pn)*x_p(pn))/rp
                  else
                    vrp=0.0d0
                    vtp=0.0d0
                  endif
! store (clock-wise) angular momentum and mass R<Rp
                  do j=idr,ndr2d-1
                    amtr(j)=amtr(j)+vtp*rp*m_p(pn)
                    mtr(j)=mtr(j)+m_p(pn)
                  enddo
                  vtzp=vz_p(pn)                 
! store data for mean and mean^2
                  vrm(idr)=vrm(idr)+vrp*m_p(pn)
                  vr2m(idr)=vr2m(idr)+(vrp**2)*m_p(pn)
                  vtm(idr)=vtm(idr)+vtp*m_p(pn)
                  vt2m(idr)=vt2m(idr)+(vtp**2)*m_p(pn)
                  vzm(idr)=vzm(idr)+vtzp*m_p(pn)
                  vz2m(idr)=vz2m(idr)+(vtzp**2)*m_p(pn)
                  z2m(idr)=z2m(idr)+(z_p(pn)**2)*m_p(pn)
                endif
              endif
            enddo
          else if(idspc(idc).eq.1) then
! for DM
            do i=0,ndm-1
              pn=i
              if(id_dm(pn).ge.ids(idc).and.id_dm(pn).le.ide(idc)) then
                rp=dsqrt(x_dm(pn)**2+y_dm(pn)**2)
                idr=int((rp-rmin2d)/dr)
                if(idr.ge.0.and.idr.lt.ndr2d) then
                  nprs=nprs+1
                  dmr(idr)=dmr(idr)+m_dm(pn)
! calculate velocity components
                  if(rp.gt.0.0d0) then
                    vrp=(vx_dm(pn)*x_dm(pn)+vy_dm(pn)*y_dm(pn))/rp
! clock-wise
                    vtp=(vx_dm(pn)*y_dm(pn)-vy_dm(pn)*x_dm(pn))/rp
                  else
                    vrp=0.0d0
                    vtp=0.0d0
                  endif
! store (clock-wise) angular momentum R<Rp
                  do j=idr,ndr2d-1
                    amtr(j)=amtr(j)+vtp*rp*m_dm(pn)
                    mtr(j)=mtr(j)+m_dm(pn)
                  enddo
                  vtzp=vz_dm(pn)                 
! store data for mean and mean^2
                  vrm(idr)=vrm(idr)+vrp*m_dm(pn)
                  vr2m(idr)=vr2m(idr)+(vrp**2)*m_dm(pn)
                  vtm(idr)=vtm(idr)+vtp*m_dm(pn)
                  vt2m(idr)=vt2m(idr)+(vtp**2)*m_dm(pn)
                  vzm(idr)=vzm(idr)+vtzp*m_dm(pn)
                  vz2m(idr)=vz2m(idr)+(vtzp**2)*m_dm(pn)
                  z2m(idr)=z2m(idr)+(z_dm(pn)**2)*m_dm(pn)
                endif
              endif
            enddo
          else
! for star
!            open(60,file='allstar.dat',status='unknown')
!            write(6,*) ' age range=',agel(idc),ageh(idc),idc
            do i=ng,ng+ns-1
              pn=list_ap(i)
!              write(60,'(4(1pE15.6))') x_p(pn),y_p(pn),z_p(pn) &
!                 ,(tu-ts_p(pn))*TMUGYR
              if(id_p(pn).ge.ids(idc).and.id_p(pn).le.ide(idc)) then
!              if(dabs(z_p(pn)).lt.0.002) then
! age constraint
              if((tu-ts_p(pn))*TMUGYR.gt.agel(idc).and. &
                (tu-ts_p(pn))*TMUGYR.lt.ageh(idc)) then
                rp=dsqrt(x_p(pn)**2+y_p(pn)**2)
                idr=int((rp-rmin2d)/dr)
                if(idr.ge.0.and.idr.lt.ndr2d) then
                  nprs=nprs+1
                  dmr(idr)=dmr(idr)+m_p(pn)
! calculate velocity components
                  if(rp.gt.0.0d0) then
                    vrp=(vx_p(pn)*x_p(pn)+vy_p(pn)*y_p(pn))/rp
! clock-wise
                    vtp=(vx_p(pn)*y_p(pn)-vy_p(pn)*x_p(pn))/rp
                  else
                    vrp=0.0d0
                    vtp=0.0d0
                  endif
! store (clock-wise) angular momentum R<Rp
                  do j=idr,ndr2d-1
                    amtr(j)=amtr(j)+vtp*rp*m_p(pn)
                    mtr(j)=mtr(j)+m_p(pn)
                  enddo
                  vtzp=vz_p(pn)                 
! store data for mean and mean^2
                  vrm(idr)=vrm(idr)+vrp*m_p(pn)
                  vr2m(idr)=vr2m(idr)+(vrp**2)*m_p(pn)
                  vtm(idr)=vtm(idr)+vtp*m_p(pn)
                  vt2m(idr)=vt2m(idr)+(vtp**2)*m_p(pn)
                  vzm(idr)=vzm(idr)+vtzp*m_p(pn)
                  vz2m(idr)=vz2m(idr)+(vtzp**2)*m_p(pn)
                  z2m(idr)=z2m(idr)+(z_p(pn)**2)*m_p(pn)
                endif
              endif
!              endif
              endif
            enddo
!            close(60)
          endif

! sum up 
! total number of particles contributed in the region of interet
            nprst=0
            call MPI_ALLREDUCE(nprs,nprst,1,MPI_INTEGER &
             ,MPI_SUM,MPI_COMM_WORLD,ierr)
            if(myrank.eq.0) then
              write(6,*) ' total number of particles within r range=',nprst
            endif

          if(nprocs.gt.1) then
            nvpr=10
            nval=nvpr*ndr2d

            allocate(tdvs(0:nval-1))            
            allocate(tdvr(0:nval-1))            

            do i=0,ndr2d-1
              tdvs(i*nvpr)=dmr(i)
              tdvs(i*nvpr+1)=vrm(i)
              tdvs(i*nvpr+2)=vr2m(i)
              tdvs(i*nvpr+3)=vtm(i)
              tdvs(i*nvpr+4)=vt2m(i)
              tdvs(i*nvpr+5)=vzm(i)
              tdvs(i*nvpr+6)=vz2m(i)
              tdvs(i*nvpr+7)=amtr(i)
              tdvs(i*nvpr+8)=mtr(i)
              tdvs(i*nvpr+9)=z2m(i)
            enddo
            do i=0,nval-1
              tdvr(i)=0.0d0
            enddo
            call MPI_ALLREDUCE(tdvs,tdvr,nval,MPI_DOUBLE_PRECISION &
             ,MPI_SUM,MPI_COMM_WORLD,ierr)
            do i=0,ndr2d-1
              dmr(i)=tdvr(i*nvpr)
              vrm(i)=tdvr(i*nvpr+1)
              vr2m(i)=tdvr(i*nvpr+2)
              vtm(i)=tdvr(i*nvpr+3)
              vt2m(i)=tdvr(i*nvpr+4)
              vzm(i)=tdvr(i*nvpr+5)
              vz2m(i)=tdvr(i*nvpr+6)
              amtr(i)=tdvr(i*nvpr+7)
              mtr(i)=tdvr(i*nvpr+8)
              z2m(i)=tdvr(i*nvpr+9)
            enddo

            deallocate(tdvr)
            deallocate(tdvs)

          endif
! calculate mean values and density
          do i=0,ndr2d-1
            if(dmr(i).gt.0.0d0) then
              vrm(i)=vrm(i)/dmr(i)
              vr2m(i)=vr2m(i)/dmr(i)
              vtm(i)=vtm(i)/dmr(i)
              vt2m(i)=vt2m(i)/dmr(i)
              vzm(i)=vzm(i)/dmr(i)
              vz2m(i)=vz2m(i)/dmr(i)
              z2m(i)=z2m(i)/dmr(i)
            endif
            ri=rmin2d+dr*dble(i)
            ro=rmin2d+dr*dble(i+1)
            dvr=M_PI*(ro**2-ri**2)
            dmr(i)=dmr(i)/dvr
          enddo

! open output file
          if(myrank.eq.0) then
            write(filename,'(a11,i1,a1,i6.6,a4)') 'output/prof',idc &
             ,'-',step,'.dat'
            open(60,file=filename,status='unknown')
            write(6,*) '# idspc,idrange=',idspc(idc),ids(idc),ide(idc)
            do i=0,ndr2d-1
              vrsig=dsqrt(vr2m(i)-vrm(i)**2)
              vtsig=dsqrt(vt2m(i)-vtm(i)**2)
              vzsig=dsqrt(vz2m(i)-vzm(i)**2)
              zhsig=dsqrt(z2m(i))
              rdmr=rmin2d+dr*(dble(i)+0.5)
              write(60,160) rdmr,dmr(i) &
               ,rdmr*LUKPC,vrsig*VUKMS,vtsig*VUKMS,vzsig*VUKMS &
               ,vtm(i)*VUKMS,zhsig*LUKPC &
               ,(rdmr+0.5d0*dr)*LUKPC,mtr(i),amtr(i)
 160          format(11(1pE13.5))
            enddo
            close(60)
          endif
        enddo
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program prof
