! /***************************************
!   bin2asc
!  3 Oct. 2017  written by D.Kawata
! ****************************************/

program bin2asc
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
      use gcdp_kernel
      use particle
      use mesh
   
      implicit none
      include 'mpif.h'

      integer i,is,nval
      integer npt,ngt,nst,ndmt,ndm1t
      integer ng,ndm,ndm1,ns
      integer nstep,step
      integer nskip,flagr,flagcom,flagbo,flagout
      integer ierr
! for checking test particle informatoin
      integer,allocatable :: flagrp(:),flagra(:),flagzm(:)
      double precision,allocatable :: rperi(:),rapo(:),zmax(:)
      double precision,allocatable :: rperi0(:),rapo0(:),zmax0(:)
! for work
      integer,allocatable :: tivr(:)
      double precision,allocatable :: tdvr(:)
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! ***  for filename ***
      character filename*60,fileo*60

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

! *****   Open Input File   ******
      if(myrank.eq.0) then
        open(50,file='./ini/input.dat',status='old')
! *****   Step number Parameter   ******
! *** fof using flag =<0: gas, 1: star, otherwise DM ***
        read(50,*) nstep,step
        read(50,*) flagr
        read(50,*) nskip
        read(50,*) flagcom
        read(50,*) flagbo
        read(50,*) flagout
        close(50)

        if(nstep.gt.1) then
          write(6,*) ' number of steps =',nstep
        else
          write(6,*) ' step =',step
        endif
        write(6,*) ' flagr = ',flagr
        write(6,*) ' nskip = ',nskip
        if(flagcom.eq.0) then
          write(6,*) ' output comoving position'
        endif
        if(flagbo.eq.0) then
          write(6,*) ' output binary for DM data (only works for DM)'
        endif
        if(flagout.eq.1) then
          write(6,*) ' output test particle info.'
        endif
      endif

      nval=7
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=nstep
        tivr(1)=step
        tivr(2)=flagr
        tivr(3)=nskip
        tivr(4)=flagcom
        tivr(5)=flagbo
        tivr(6)=flagout
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nstep=tivr(0)
      step=tivr(1)
      flagr=tivr(2)
      nskip=tivr(3)
      flagcom=tivr(4)
      flagbo=tivr(5)
      flagout=tivr(6)
      deallocate(tivr)

      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      if(nstep.gt.1.and.myrank.eq.0) then
        open(49,file='./ini/file.dat',status='old')
      endif
      do is=0,nstep-1
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

        if(is.eq.0) then
          allocate(flagrp(0:ndm-1))
          allocate(flagra(0:ndm-1))
          allocate(flagzm(0:ndm-1))
          allocate(rperi(0:ndm-1))
          allocate(rapo(0:ndm-1))
          allocate(zmax(0:ndm-1))
          allocate(rperi0(0:ndm-1))
          allocate(rapo0(0:ndm-1))
          allocate(zmax0(0:ndm-1))
        endif

        if(flagout.eq.1) then
          if(is.eq.0) then
            do i=0,ndm-1
              flagrp(i)=0
              flagra(i)=0
              flagzm(i)=0
              rperi(i)=rperi_dm(i)
              rapo(i)=rapo_dm(i)
              zmax(i)=zmax_dm(i)
              rperi0(i)=rperi_dm(i)
              rapo0(i)=rapo_dm(i)
              zmax0(i)=zmax_dm(i)
            enddo
          else
            do i=0,ndm-1
              if(flagrp(i).eq.0) then
                if(rperi(i).eq.rperi_dm(i).and.rperi0(i).ne.rperi_dm(i)) then 
                  flagrp(i)=1
                else
                  rperi(i)=rperi_dm(i)
                endif
              else
                rperi_dm(i)=rperi(i)
              endif
              if(flagra(i).eq.0) then
                if(rapo(i).eq.rapo_dm(i).and.rapo0(i).ne.rapo_dm(i)) then 
                  flagra(i)=1
                else
                  rapo(i)=rapo_dm(i)
                endif
              else
                rapo_dm(i)=rapo(i)
              endif
              if(flagzm(i).eq.0) then
                if(zmax(i).eq.zmax_dm(i).and.zmax0(i).ne.zmax_dm(i)) then 
                  flagzm(i)=1
                else
                  zmax(i)=zmax_dm(i)
                endif
              else
                zmax_dm(i)=zmax(i)
              endif
            enddo
          endif
        endif

        if(flagcom.ne.0) then
! output physical value
          ai=1.0d0
        endif
        call output(ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns &
          ,flagr,step,nskip,ai,flagbo,flagout)

      enddo
      if(myrank.eq.0) then
        close(49)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program
