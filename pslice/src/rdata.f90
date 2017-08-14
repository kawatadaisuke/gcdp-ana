! *************************************************************
!    rdatav.f  read data for the data from the code > ver.32.3
!  14 Jan. 2016  written by D.Kawata
! **************************************************************
!70*******************************************************************
! flago!=0 case is not done yet.
!

subroutine rdata(step,ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns &
       ,ai,tn,flagr)
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm

      include 'mpif.h'

      integer i,j,flagr
      integer step,ngt,ng,ndmt,ndm,nst,ns,ndm1t,ndm1,npt,np,npst
      integer ngr,ndmr,nsr,ndm1r,nagr,nadmr,nasr
      integer npr,nc,flagf,isend,pn
      integer nof,ifn,ip,cndmt,cngt,cnst,iv
      integer ndval,nival,nprocr,ndbhyd,ndbmet,ndbsf
      double precision ai,tn
      integer ng0,ns0
      character filei*60
! *** for work ***
      integer invali,invald,nval
      integer,allocatable :: idisp(:),jjlen(:)
      integer,allocatable :: npstp(:),npenp(:),jstas(:),jends(:)
      integer,allocatable :: tivr(:),tivs(:),lists(:)
      double precision,allocatable :: tdvr(:),tdvs(:)

      allocate(jstas(0:nprocs-1))
      allocate(jends(0:nprocs-1))
      allocate(idisp(0:nprocs-1))
      allocate(jjlen(0:nprocs-1))

      write(6,*) ' step,flagr=',step,flagr

! number of particle for each proc
      ng=0
      ndm=0
      ndm1=0
      ns=0
      np=0
      npt=0
      ngt=0
      nst=0

      ifn=0
      npst=0
      cndmt=0
      cngt=0
      cnst=0
 70   flagf=0
      if(myrank.eq.0) then
        write(filei,'(a21,i6.6,a1,i4.4)') '../output/data/bdvals',step,'n',ifn
        write(6,*) ' reading ',filei
        open(50,file=filei,status='old',form='unformatted',err=90)
        read(50) npt,ndmt,ndm1t,ai,tn
        read(50) nprocr,nof,invali,invald
        if(ifn.eq.0) then
          write(6,*) ' npt,ndmt,ndm1t,nprocr,nof=',npt,ndmt,ndm1t,nprocr,nof
          write(6,*) ' nvali,nvald=',invali,invald
        endif
        goto 91
 90     flagf=1
      endif
 91   call MPI_BCAST(flagf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(flagf.ne.0) then
! no DM data
        goto 92
      endif

      nval=6
      allocate(tivr(0:nval-1))

      if(myrank.eq.0) then
        tivr(0)=ndmt
        tivr(1)=ndm1t
        tivr(2)=nprocr
        tivr(3)=nof
        tivr(4)=invali
        tivr(5)=invald
      endif
      call MPI_BCAST(tivr,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        ndmt=tivr(0)
        ndm1t=tivr(1)
        nprocr=tivr(2)
        nof=tivr(3)
        invali=tivr(4)
        invald=tivr(5)
      endif

      deallocate(tivr)

! nprocr: number of proc used in simualtion. 
      if(allocated(npstp)) then
        deallocate(npstp)
        deallocate(npenp)
      endif
      allocate(npstp(0:nprocr-1))
      allocate(npenp(0:nprocr-1))

! number of particles for each proc
      if(myrank.eq.0) then
        do i=0,nprocr-1
          npstp(i)=0
          npenp(i)=0
        enddo
      endif
! *** get start and end id particles for each core recieve ***
      do i=0,nprocs-1
        call para_range(0,ndmt-1,nprocs,i,jsta,jend)  
        jstas(i)=jsta
        jends(i)=jend
      enddo

      if(myrank.eq.0) then
        do ip=ifn,nprocr-1,nof
          read(50) ngr,ndmr,nsr,ndm1r,nagr,nadmr,nasr
!          write(6,*) ifn,ip,' ng,ndm,ns=',ngr,ndmr,nsr,cndmt
          cndmt=cndmt+ndmr
          npstp(ip)=npst
          npenp(ip)=npst+ndmr-1
          npst=npenp(ip)+1
        enddo
      endif
      call MPI_BCAST(npstp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npenp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      nc=ndm
      do ip=ifn,nprocr-1,nof
! *** set idisp and jjlen ***
        do i=0,nprocs-1
          idisp(i)=0
          jjlen(i)=0
        enddo
! *** count number of particles read for each proc
        do i=npstp(ip),npenp(ip)
          do j=0,nprocs-1
            if(i.ge.jstas(j).and.i.le.jends(j)) then
              if(jjlen(j).eq.0) then
                idisp(j)=isend
              endif
              jjlen(j)=jjlen(j)+1
            endif
          enddo
        enddo
        nc=nc+jjlen(myrank)
      enddo
      ndm=ndm+nc
      ifn=ifn+1
!      if(myrank.eq.0) then
!        write(6,*) ifn,'file, cndmt,ndm,nc,jjlen=',cndmt,ndm,nc,jjlen(myrank)
!      endif
      close(50)
      if(ifn.lt.nof) then
        goto 70
      endif

      call allocate_dm(ndm)
!      if(myrank.eq.0) then
!        write(6,*) ' for counting: cndmt,ndm=',cndmt,ndm,nc,jjlen(myrank)
!      endif

! reset ndm and etc.
      ndm=0
      cndmt=0
      npst=0
      ifn=0


 72   if(myrank.eq.0) then
        write(filei,'(a21,i6.6,a1,i4.4)') '../output/data/bdvals',step,'n',ifn
!        write(6,*) ' reading ',filei
        open(50,file=filei,status='old',form='unformatted')
        read(50) npt,ndmt,ndm1t,ai,tn
        read(50) nprocr,nof,invali,invald
      endif

      nval=6
      allocate(tivr(0:nval-1))

      if(myrank.eq.0) then
        tivr(0)=ndmt
        tivr(1)=ndm1t
        tivr(2)=nprocr
        tivr(3)=nof
        tivr(4)=invali
        tivr(5)=invald
      endif
      call MPI_BCAST(tivr,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        ndmt=tivr(0)
        ndm1t=tivr(1)
        nprocr=tivr(2)
        nof=tivr(3)
        invali=tivr(4)
        invald=tivr(5)
      endif

      deallocate(tivr)

! nprocr: number of proc used in simualtion. 
      if(allocated(npstp)) then
        deallocate(npstp)
        deallocate(npenp)
      endif
      allocate(npstp(0:nprocr-1))
      allocate(npenp(0:nprocr-1))

! number of particles for each proc
      if(myrank.eq.0) then
        do i=0,nprocr-1
          npstp(i)=0
          npenp(i)=0
        enddo
      endif
! *** get start and end id particles for each core recieve ***
      do i=0,nprocs-1
        call para_range(0,ndmt-1,nprocs,i,jsta,jend)  
        jstas(i)=jsta
        jends(i)=jend
      enddo

      if(myrank.eq.0) then
        do ip=ifn,nprocr-1,nof
          read(50) ngr,ndmr,nsr,ndm1r,nagr,nadmr,nasr
!          write(6,*) ifn,ip,' ng,ndm,ns=',ngr,ndmr,nsr,cndmt
          cndmt=cndmt+ndmr
          npstp(ip)=npst
          npenp(ip)=npst+ndmr-1
          npst=npenp(ip)+1
        enddo
      endif
      call MPI_BCAST(npstp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npenp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! *** reading integer data ***
      do iv=1,invali
        nc=ndm
        do ip=ifn,nprocr-1,nof
          npr=npenp(ip)-npstp(ip)+1

          allocate(tivs(0:npr-1))

          if(myrank.eq.0) then
            read(50) (tivs(i),i=0,npr-1)
          endif
! *** set idisp and jjlen ***
          do i=0,nprocs-1
            idisp(i)=0
            jjlen(i)=0
          enddo
! *** count number of particles read for each proc
          isend=0
          do i=npstp(ip),npenp(ip)
            do j=0,nprocs-1
              if(i.ge.jstas(j).and.i.le.jends(j)) then
                if(jjlen(j).eq.0) then
                  idisp(j)=isend
                endif
                jjlen(j)=jjlen(j)+1
              endif
            enddo
            isend=isend+1
          enddo
          npr=jjlen(myrank)

          allocate(tivr(0:npr-1))

          call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER &
            ,tivr,npr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          if(iv.eq.1) then
            do i=0,jjlen(myrank)-1
              id_dm(i+nc)=tivr(i)
            enddo
          else if(iv.eq.2) then
            do i=0,jjlen(myrank)-1
              list_adm(i+nc)=tivr(i)
            enddo
          endif
          nc=nc+jjlen(myrank)

          deallocate(tivs)
          deallocate(tivr)

        enddo
      enddo

! *** reading double precision data ***
      do iv=1,invald
        nc=ndm
        do ip=ifn,nprocr-1,nof
          npr=npenp(ip)-npstp(ip)+1

          allocate(tdvs(0:npr-1))

          if(myrank.eq.0) then
            read(50) (tdvs(i),i=0,npr-1)
          endif
! *** set idisp and jjlen ***
          do i=0,nprocs-1
            idisp(i)=0
            jjlen(i)=0
          enddo
! *** count number of particles read for each proc
          isend=0
          do i=npstp(ip),npenp(ip)
            do j=0,nprocs-1
              if(i.ge.jstas(j).and.i.le.jends(j)) then
                if(jjlen(j).eq.0) then
                  idisp(j)=isend
                endif
                jjlen(j)=jjlen(j)+1
              endif
            enddo
            isend=isend+1
          enddo
!          write(6,*) ' idisp,jjlen=',idisp(0),jjlen(0) &
!            ,jstas(0),jends(0),npstp(ip),npenp(ip),ip
          npr=jjlen(myrank)

          allocate(tdvr(0:npr-1))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,tdvr,npr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          if(iv.eq.1) then
            do i=0,jjlen(myrank)-1
              x_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.2) then
            do i=0,jjlen(myrank)-1
              y_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.3) then
            do i=0,jjlen(myrank)-1
              z_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.4) then
            do i=0,jjlen(myrank)-1
              vx_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.5) then
            do i=0,jjlen(myrank)-1
              vy_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.6) then
            do i=0,jjlen(myrank)-1
              vz_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.7) then
            do i=0,jjlen(myrank)-1
              m_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.8) then
            do i=0,jjlen(myrank)-1
              rho_dm(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.9) then
            do i=0,jjlen(myrank)-1
              h_dm(i+nc)=tdvr(i)
            enddo
          endif
          nc=nc+jjlen(myrank)

          deallocate(tdvs)
          deallocate(tdvr)

        enddo
      enddo
      ndm=ndm+nc
      ifn=ifn+1
!      if(myrank.eq.0) then
!        write(6,*) ' cndmt,ndm=',cndmt,ndm,nc,jjlen(myrank)
!      endif
      close(50)
      if(ifn.lt.nof) then
        goto 72
      endif
      if(myrank.eq.0.and.cndmt.ne.ndmt) then
        write(6,*) ' Number of total DM particle is inconsistent'
        write(6,*) ' ndmt (read)=',cndmt
        write(6,*) ' ndmt (file header) =',ndmt
        call MPI_FINALIZE()
        stop
      endif
      if(ndm.ne.jends(myrank)-jstas(myrank)+1) then
        if(myrank.eq.0) then
          write(6,*) ' NDM particle in each proc is inconsistent'
          write(6,*) ' ndm (read)=',ndm
          write(6,*) ' ndm (from npt)=',jends(myrank)-jstas(myrank)+1
        endif
        call MPI_FINALIZE()
        stop
      endif

! *** read baryon data ***
 92   ifn=0
      ngt=0
      npst=0
      np=0

! count baryon particle in each core
 71   flagf=0
      if(myrank.eq.0) then
        write(filei,'(a21,i6.6,a1,i4.4)') &
         '../output/data/bbvals',step,'n',ifn
        open(50,file=filei,status='old',form='unformatted',err=93)
        read(50) npt,ndmt,ndm1t,ai,tn
        read(50) nprocr,nof,invali,invald
        if(ifn.eq.0) then
          write(6,*) ' bbvals*: nvali,nvald=',invali,invald
        endif
        goto 94
 93     flagf=1
      endif
 94   call MPI_BCAST(flagf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(flagf.ne.0) then
! no baryon data
        goto 95
      endif

      allocate(tivr(0:5))

      if(myrank.eq.0) then
        tivr(0)=npt
        tivr(1)=ndmt
        tivr(2)=nprocr
        tivr(3)=nof
        tivr(4)=invali
        tivr(5)=invald
      endif
      call MPI_BCAST(tivr,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        npt=tivr(0)
        ndmt=tivr(1)
        nprocr=tivr(2)
        nof=tivr(3)
        invali=tivr(4)
        invald=tivr(5)
      endif

      deallocate(tivr)
! nprocr: number of proc used in simualtion. 
      if(allocated(npstp)) then
        deallocate(npstp)
        deallocate(npenp)
      endif
      allocate(npstp(0:nprocr-1))
      allocate(npenp(0:nprocr-1))

! number of particles for each proc
      if(myrank.eq.0) then
        do i=0,nprocr-1
          npstp(i)=0
          npenp(i)=0
        enddo
      endif
! *** get start and end id particles for each core recieve ***
      do i=0,nprocs-1
        call para_range(0,npt-1,nprocs,i,jsta,jend)  
        jstas(i)=jsta
        jends(i)=jend
      enddo

      if(myrank.eq.0) then
        do ip=ifn,nprocr-1,nof
           read(50) ngr,ndmr,nsr,ndm1r,nagr,nadmr,nasr
!        write(6,*) ifn,ip,' ng,ndm,ns=',ng,ndm,ns
          npstp(ip)=npst
          npenp(ip)=npst+ngr+nsr-1
          npst=npenp(ip)+1        
        enddo
      endif
      call MPI_BCAST(npstp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npenp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nc=np
      do ip=ifn,nprocr-1,nof
! *** set idisp and jjlen ***
        do i=0,nprocs-1
          idisp(i)=0
          jjlen(i)=0
        enddo
! *** count number of particles read for each proc
        do i=npstp(ip),npenp(ip)
          do j=0,nprocs-1
            if(i.ge.jstas(j).and.i.le.jends(j)) then
              if(jjlen(j).eq.0) then
                idisp(j)=isend
              endif
              jjlen(j)=jjlen(j)+1
            endif
          enddo
          isend=isend+1
        enddo
        nc=nc+jjlen(myrank)
      enddo
      np=np+nc
      ifn=ifn+1
      close(50)
      if(ifn.lt.nof) then
        goto 71
      endif

      call allocate_baryon(np)
!      if(myrank.eq.0) then
!        write(6,*) ' for counting np(rank=0)=',np,myrank
!      endif

      ifn=0
      ngt=0
      npst=0
      np=0

 73   if(myrank.eq.0) then
        write(filei,'(a21,i6.6,a1,i4.4)') &
         '../output/data/bbvals',step,'n',ifn
        open(50,file=filei,status='old',form='unformatted')
        read(50) npt,ndmt,ndm1t,ai,tn
        read(50) nprocr,nof,invali,invald
        if(ifn.eq.0) then
          write(6,*) ' bbvals*: nvali,nvald=',invali,invald
        endif
      endif
      if(flagr.gt.0) then
! reading extra data
        write(filei,'(a21,i6.6,a1,i4.4)') &
        '../output/data/bbhyds',step,'n',ifn
        open(51,file=filei,status='old',form='unformatted',err=96)
        read(51) npt,ndmt,ndm1t,ai,tn
        read(51) nprocr,nof,ndbhyd

        write(filei,'(a21,i6.6,a1,i4.4)') &
          '../output/data/bbmets',step,'n',ifn
        open(52,file=filei,status='old',form='unformatted',err=96)
        read(52) npt,ndmt,ndm1t,ai,tn
        read(52) nprocr,nof,ndbmet

        write(filei,'(a21,i6.6,a1,i4.4)') &
        '../output/data/bbsfis',step,'n',ifn
        open(54,file=filei,status='old',form='unformatted',err=96)
        read(54) npt,ndmt,ndm1t,ai,tn
        read(54) nprocr,nof,ndbsf

      else
        ndbhyd=0
        ndbmet=0
        nbdsf=0
      endif

 96   allocate(tivr(0:6))

      if(myrank.eq.0) then
        tivr(0)=npt
        tivr(1)=ndmt
        tivr(2)=nprocr
        tivr(3)=nof
        tivr(4)=invali
        tivr(5)=invald
        tivr(6)=ndbhyd
      endif
      call MPI_BCAST(tivr,7,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myrank.ne.0) then
        npt=tivr(0)
        ndmt=tivr(1)
        nprocr=tivr(2)
        nof=tivr(3)
        invali=tivr(4)
        invald=tivr(5)
        ndbhyd=tivr(6)
      endif

      deallocate(tivr)
! nprocr: number of proc used in simualtion. 
      if(allocated(npstp)) then
        deallocate(npstp)
        deallocate(npenp)
      endif
      allocate(npstp(0:nprocr-1))
      allocate(npenp(0:nprocr-1))

! number of particles for each proc
      if(myrank.eq.0) then
        do i=0,nprocr-1
          npstp(i)=0
          npenp(i)=0
        enddo
      endif
! *** get start and end id particles for each core recieve ***
      do i=0,nprocs-1
        call para_range(0,npt-1,nprocs,i,jsta,jend)  
        jstas(i)=jsta
        jends(i)=jend
      enddo

      if(myrank.eq.0) then
        do ip=ifn,nprocr-1,nof
           read(50) ngr,ndmr,nsr,ndm1r,nagr,nadmr,nasr
          npstp(ip)=npst
          npenp(ip)=npst+ngr+nsr-1
          npst=npenp(ip)+1        
        enddo
      endif
      call MPI_BCAST(npstp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npenp,nprocr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! *** initialization ***
      do i=0,np-1
        flagfd_p(i)=0
      enddo
      do iv=1,invali
        nc=np
        do ip=ifn,nprocr-1,nof
          npr=npenp(ip)-npstp(ip)+1

          allocate(tivs(0:npr-1))

          if(myrank.eq.0) then
            read(50) (tivs(i),i=0,npr-1)
          endif
! *** set idisp and jjlen ***
          do i=0,nprocs-1
            idisp(i)=0
            jjlen(i)=0
          enddo
! *** count number of particles read for each proc
          isend=0
          do i=npstp(ip),npenp(ip)
            do j=0,nprocs-1
              if(i.ge.jstas(j).and.i.le.jends(j)) then
                if(jjlen(j).eq.0) then
                  idisp(j)=isend
                endif
                jjlen(j)=jjlen(j)+1
              endif
            enddo
            isend=isend+1
          enddo
          npr=jjlen(myrank)

          allocate(tivr(0:npr-1))

          call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER &
           ,tivr,npr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          if(iv.eq.1) then
            do i=0,jjlen(myrank)-1
              id_p(i+nc)=tivr(i)
            enddo
          else if(iv.eq.2) then
            do i=0,jjlen(myrank)-1
              flagc_p(i+nc)=tivr(i)
            enddo
          else if(iv.eq.3) then
            do i=0,jjlen(myrank)-1
              list_ap(i+nc)=tivr(i)
            enddo
          else if(iv.eq.4) then
            do i=0,jjlen(myrank)-1
              flagfd_p(i+nc)=tivr(i)
            enddo
          endif
          nc=nc+jjlen(myrank)

          deallocate(tivs)
          deallocate(tivr)

        enddo
      enddo
! *** reading double precision data
      do iv=1,invald
        nc=np
        do ip=ifn,nprocr-1,nof
          npr=npenp(ip)-npstp(ip)+1

          allocate(tdvs(0:npr-1))

          if(myrank.eq.0) then
            read(50) (tdvs(i),i=0,npr-1)
          endif
! *** set idisp and jjlen ***
          do i=0,nprocs-1
            idisp(i)=0
            jjlen(i)=0
          enddo
! *** count number of particles read for each proc
          isend=0
          do i=npstp(ip),npenp(ip)
            do j=0,nprocs-1
              if(i.ge.jstas(j).and.i.le.jends(j)) then
                if(jjlen(j).eq.0) then
                  idisp(j)=isend
                endif
                jjlen(j)=jjlen(j)+1
              endif
            enddo
            isend=isend+1
          enddo
          npr=jjlen(myrank)

          allocate(tdvr(0:npr-1))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,tdvr,npr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          if(iv.eq.1) then
            do i=0,jjlen(myrank)-1
              x_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.2) then
            do i=0,jjlen(myrank)-1
              y_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.3) then
            do i=0,jjlen(myrank)-1
              z_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.4) then
            do i=0,jjlen(myrank)-1
              vx_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.5) then
            do i=0,jjlen(myrank)-1
              vy_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.6) then
            do i=0,jjlen(myrank)-1
              vz_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.7) then
            do i=0,jjlen(myrank)-1
              m_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.8) then
            do i=0,jjlen(myrank)-1
              rho_p(i+nc)=tdvr(i)
            enddo
          else if(iv.eq.9) then
            do i=0,jjlen(myrank)-1
              u_p(i+nc)=tdvr(i)
            enddo
          endif
          nc=nc+jjlen(myrank)

          deallocate(tdvs)
          deallocate(tdvr)

        enddo
      enddo
      if(flagr.gt.0) then
! *** hyd vals
        do iv=1,ndbhyd
          nc=np
          do ip=ifn,nprocr-1,nof
            npr=npenp(ip)-npstp(ip)+1

            allocate(tdvs(0:npr-1))

            if(myrank.eq.0) then
              read(51) (tdvs(i),i=0,npr-1)
            endif
! *** set idisp and jjlen ***
            do i=0,nprocs-1
              idisp(i)=0
              jjlen(i)=0
            enddo
! *** count number of particles read for each proc
            isend=0
            do i=npstp(ip),npenp(ip)
              do j=0,nprocs-1
                if(i.ge.jstas(j).and.i.le.jends(j)) then
                  if(jjlen(j).eq.0) then
                    idisp(j)=isend
                  endif
                  jjlen(j)=jjlen(j)+1
                endif
              enddo
              isend=isend+1
            enddo
            npr=jjlen(myrank)

            allocate(tdvr(0:npr-1))

            call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
             ,tdvr,npr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            if(iv.eq.1) then
              do i=0,jjlen(myrank)-1
                h_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.2) then
              do i=0,jjlen(myrank)-1
                div_v_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.3) then
              do i=0,jjlen(myrank)-1
                alpv_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.4) then
              do i=0,jjlen(myrank)-1
                alpu_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.5) then
              do i=0,jjlen(myrank)-1
                myu_p(i+nc)=tdvr(i)
              enddo
            endif
            nc=nc+jjlen(myrank)

            deallocate(tdvs)
            deallocate(tdvr)

          enddo
        enddo
! metal values
        do iv=1,ndbmet
          nc=np
          do ip=ifn,nprocr-1,nof
            npr=npenp(ip)-npstp(ip)+1

            allocate(tdvs(0:npr-1))

            if(myrank.eq.0) then
              read(52) (tdvs(i),i=0,npr-1)
            endif
! *** set idisp and jjlen ***
            do i=0,nprocs-1
              idisp(i)=0
              jjlen(i)=0
            enddo
! *** count number of particles read for each proc
            isend=0
            do i=npstp(ip),npenp(ip)
              do j=0,nprocs-1
                if(i.ge.jstas(j).and.i.le.jends(j)) then
                  if(jjlen(j).eq.0) then
                    idisp(j)=isend
                  endif
                  jjlen(j)=jjlen(j)+1
                endif
              enddo
              isend=isend+1
            enddo
            npr=jjlen(myrank)

            allocate(tdvr(0:npr-1))

            call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
             ,tdvr,npr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            if(iv.eq.1) then
              do i=0,jjlen(myrank)-1
                mzHe_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.2) then
              do i=0,jjlen(myrank)-1
                mzZ_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.3) then
              do i=0,jjlen(myrank)-1
                mzC_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.4) then
              do i=0,jjlen(myrank)-1
                mzN_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.5) then
              do i=0,jjlen(myrank)-1
                mzO_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.6) then
              do i=0,jjlen(myrank)-1
                mzNe_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.7) then
              do i=0,jjlen(myrank)-1
                mzMg_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.8) then
              do i=0,jjlen(myrank)-1
                mzSi_p(i+nc)=tdvr(i)
              enddo
            else if(iv.eq.9) then
              do i=0,jjlen(myrank)-1
                mzFe_p(i+nc)=tdvr(i)
              enddo
            endif
            nc=nc+jjlen(myrank)

            deallocate(tdvs)
            deallocate(tdvr)

          enddo
        enddo
! star data
        do iv=1,ndbsf
          nc=np
          do ip=ifn,nprocr-1,nof
            npr=npenp(ip)-npstp(ip)+1

            allocate(tdvs(0:npr-1))

            if(myrank.eq.0) then
              read(54) (tdvs(i),i=0,npr-1)
            endif
! *** set idisp and jjlen ***
            do i=0,nprocs-1
              idisp(i)=0
              jjlen(i)=0
            enddo
! *** count number of particles read for each proc
            isend=0
            do i=npstp(ip),npenp(ip)
              do j=0,nprocs-1
                if(i.ge.jstas(j).and.i.le.jends(j)) then
                  if(jjlen(j).eq.0) then
                    idisp(j)=isend
                  endif
                  jjlen(j)=jjlen(j)+1
                endif
              enddo
              isend=isend+1
            enddo
            npr=jjlen(myrank)

            allocate(tdvr(0:npr-1))

            call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
             ,tdvr,npr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            if(iv.eq.1) then
              do i=0,jjlen(myrank)-1
                ts_p(i+nc)=tdvr(i)
              enddo
            endif
            nc=nc+jjlen(myrank)

            deallocate(tdvs)
            deallocate(tdvr)

          enddo
        enddo

        close(51)
        close(52)
        close(54)
      endif
! update np and ifn
      np=np+nc
      ifn=ifn+1
      close(50)
      if(ifn.lt.nof) then
        goto 73
      endif
      if(np.ne.jends(myrank)-jstas(myrank)+1) then
        if(myrank.eq.0) then
          write(6,*) ' Nbaryon particle in each proc is inconsistent',myrank
          write(6,*) ' np (read)=',np,myrank
          write(6,*) ' np (from npt)=',jends(myrank)-jstas(myrank)+1,myrank
        endif
        call MPI_FINALIZE()
        stop
      endif

! *** re-construct list_ap ***
 95   cnst=0     
      cngt=0    

      allocate(lists(0:np-1))
 
      do i=0,np-1
        if(flagc_p(i).le.0) then
          list_ap(cngt)=i
          cngt=cngt+1
        else
          lists(cnst)=i
          cnst=cnst+1
        endif
      enddo
      ng=cngt
      ns=cnst
      do i=0,ns-1
        list_ap(ng+i)=lists(i)
      enddo

      deallocate(lists)

! *** get total number of star and gas particle

      allocate(tivs(0:1))
      allocate(tivr(0:1))

      tivs(0)=ng
      tivs(1)=ns
      call MPI_ALLREDUCE(tivs,tivr,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      ngt=tivr(0)
      nst=tivr(1)

      deallocate(tivs)
      deallocate(tivr)

      if(myrank.eq.0) then
        write(6,*) ' ngt,nst,npt=',ngt,nst,npt
      endif
!      write(6,*) ' myrank,ng,ns,np=',myrank,ng,ns,np

! *** re-construct list_adm ***

     allocate(lists(0:ndm-1))

      cnst=0     
      cngt=0     
      do i=0,ndm-1
        if(id_dm(i).lt.ndm1t) then
          list_adm(cngt)=i
          cngt=cngt+1
        else
          lists(cnst)=i
          cnst=cnst+1
        endif
      enddo
      ndm1=cngt
      do i=0,ndm-ndm1-1
        list_adm(ndm1+i)=lists(i)
      enddo
      
!      write(6,*) ' myrank,ndm,ndm1=',myrank,ndm,ndm1,ndm1t

      deallocate(lists)

! *** send ai,tn

      allocate(tdvr(0:1))

      tdvr(0)=ai
      tdvr(1)=tn
      call MPI_BCAST(tdvr,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      ai=tdvr(0)
      tn=tdvr(1)

      deallocate(tdvr)

!      write(filei,'(a3,i3.3)') 'hdm',myrank
!      open(60,file=filei,status='unknown')
!      do i=0,ndm1-1
!        pn=list_adm(i)
!        write(60,'(4(1pE13.5),I10)') x_dm(pn),y_dm(pn),z_dm(pn),m_dm(pn),id_dm(pn)
!      enddo
!      close(60)

!      write(filei,'(a3,i3.3)') 'ldm',myrank
!      open(60,file=filei,status='unknown')
!      do i=ndm1,ndm-1
!        pn=list_adm(i)
!        write(60,'(4(1pE13.5),I10)') x_dm(pn),y_dm(pn),z_dm(pn),m_dm(pn),id_dm(pn)
!      enddo

!      close(60)

!      write(filei,'(a3,i3.3)') 'gas',myrank
!      open(60,file=filei,status='unknown')
!      do i=0,ng-1
!        pn=list_ap(i)
!        write(60,'(3(1pE13.5),I10)') x_p(pn),y_p(pn),z_p(pn),id_p(pn)
!      enddo
!      close(60)

!      write(filei,'(a4,i3.3)') 'star',myrank
!      open(60,file=filei,status='unknown')
!      do i=0,ns-1
!        pn=list_ap(i+ng)
!        write(60,'(3(1pE13.5),I10)') x_p(pn),y_p(pn),z_p(pn),id_p(pn)
!      enddo
!      close(60)

      deallocate(npstp)
      deallocate(npenp)
      deallocate(jstas)
      deallocate(jends)
      deallocate(idisp)
      deallocate(jjlen)


end subroutine
