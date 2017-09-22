! ***********************************************
!    ddecb.f90  for epotgrid
!  04  Dec. 2014    written by D.KAWATA
! ***********************************************

! ***********************************************
!  Domain decomposition for baryon particles
!  using Paeno-Hilbert ordering 
! ***********************************************

subroutine ddecb(npt,np,ng,ns)
      use gcdp_const
      use gcdp_baryon
      use gcdp_btree
      use gcdp_system

      implicit none
      include 'mpif.h'

      integer,intent(in) :: npt
      integer,intent(inout) :: np,ng,ns
      double precision kmmp
      integer i,j,k,level,maxlevel,level_mergin,level_limit,ip,nin,ierr
! *** diameter of subcell ***	  
      double precision l_sc
!  center of the simulation box 	  
      double precision cx0,cy0,cz0
! ** for parallel process ***
      integer npj,crank,isend
      integer nvali,nval,jsta,jend
      integer ireqs(0:3),ireqr(0:3)
      double precision dx(0:1)
      data dx/1.0d0,-1.0d0/
! *** for domain decomposition ver. 2 ***
! *** to find active particles ***
      integer pn,na,nanp
      double precision lwt,upt
! *** for check ***
      integer nptc
! *** for test ***
      character fileo*60
! for work
      integer,allocatable :: istatus(:),istatusr(:)
      integer,allocatable :: srank(:),rrank(:)
      integer,allocatable :: jstap(:),jendp(:)
      integer,allocatable :: npsp(:),nprp(:)
      integer,allocatable :: idisp(:),jjlen(:)
      integer,allocatable :: ord(:),ord1(:),pn_nfp(:),nd_nfp(:)
      integer,allocatable :: list(:),talist(:),nalist(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)
! ** for Paeno-Hilbert ordering **
      integer,allocatable :: ixp(:),iyp(:),izp(:),snfp(:),c_nfp(:)
      integer*8,allocatable :: phkey(:)
      integer*8,allocatable :: phkeyr(:),phkeys(:)
! subcell
      integer numsc,pnumsc,nscnf
      integer,allocatable :: npsc(:),dpsc(:),scp(:),nppsc(:)
      integer*8,allocatable :: phkey1(:)

! *** mergin of tree level to make a key for the particles ***
      level_mergin = 1
! *** phkey < 2**64 ***
      level_limit=20-level_mergin

!       write(fileo,'(a5,i3.3)') 'pdecb',myrank
!       open(60,file=fileo,status='unknown')
!       write(60,*) '# np=',np
!      do i=0,np-1
!       write(60,'(5(1pE13.5),I10)') x_p(i),y_p(i),z_p(i) &
!        ,m_p(i),h_p(i),id_p(i)
!      enddo
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      close(60)

! ** Make root in tree.f **	
! *** since ver.33. need to list pn_nfp

      allocate(pn_nfp(0:np))      

      if(.not.allocated(np_tr)) then
! allocate memory for tree, need only the np size of _tr array
        call allocate_btree(np)
      endif

      do i=0,np-1
        pn_nfp(i)=i
      enddo
      call makeroot(np,pn_nfp,0.0d0)
      cx0=cx_tr(0)
      cy0=cy_tr(0)
      cz0=cz_tr(0)
      maxlevel=0

      deallocate(pn_nfp)
      allocate(phkey(0:np))
      allocate(c_nfp(0:np))
      allocate(snfp(0:np))

!  initialization 
      do i=0,np-1
        phkey(i)=0
! *** use cx_tr for the position of the cell for the particles 
        cx_tr(i)=cx0
        cy_tr(i)=cy0
        cz_tr(i)=cz0
        c_nfp(i)=0
        snfp(i)=0
      enddo
      l_sc=l_tr(0)*0.5d0	
      level = 0

! allocate daughter subcell
      numsc=NDTREE

      allocate(npsc(0:numsc-1))
      allocate(dpsc(0:numsc-1))

      do i=0,NDTREE-1
        npsc(i)=0
        dpsc(i)=0
      enddo

! record current subcell for particles
        
      allocate(scp(0:np-1))

      do i=0,np-1
        scp(i)=0
      enddo

! ****   start iteration in level ****
   77 level=level+1	  
      l_sc = 0.5d0*l_sc
!  find out which subcell it is in 

      allocate(ixp(0:np-1))
      allocate(iyp(0:np-1))
      allocate(izp(0:np-1))

      do i=0,np-1
        if(x_p(i)-cx_tr(i).ge.0.0d0) then
          ixp(i)=0 
        else
          ixp(i)=1
        endif 
        if(y_p(i)-cy_tr(i).ge.0.0d0) then
          iyp(i)=0 
        else
          iyp(i)=1
        endif 
        if(z_p(i)-cz_tr(i).ge.0.0d0) then
          izp(i)=0 
        else
          izp(i)=1
        endif 
      enddo
! *** set c and state ***
      call phcurven(0,np-1,snfp,c_nfp,ixp,iyp,izp,level)

      do i=0,np-1
! *** asign the ph key ***
! *** has to be phkey(i)*NDPH 
        phkey(i)=phkey(i)*NDPH+c_nfp(i)
! *** update cx for the cell ***
        cx_tr(i)=dx(ixp(i))*l_sc+cx_tr(i)
        cy_tr(i)=dx(iyp(i))*l_sc+cy_tr(i)
        cz_tr(i)=dx(izp(i))*l_sc+cz_tr(i)
      enddo

      deallocate(ixp)
      deallocate(iyp)
      deallocate(izp)

      do i=0,np-1
! allocate subcell
        scp(i)=dpsc(scp(i))+c_nfp(i)
! counting number of particles in subcell
        npsc(scp(i))=npsc(scp(i))+1
      enddo

! make new subcell
      pnumsc=numsc
      numsc=0
      nscnf=0
! counting numsc for the next level
      do i=0,pnumsc-1
        if(npsc(i).gt.NDTREE) then
          nscnf=nscnf+1
        endif
        if(npsc(i).gt.0) then
          do j=0,NDTREE-1
            numsc=numsc+1
          enddo
        endif
      enddo

      allocate(nppsc(0:pnumsc-1))

      do i=0,pnumsc-1
        nppsc(i)=npsc(i)
      enddo

      deallocate(dpsc)
      deallocate(npsc)

      allocate(dpsc(0:pnumsc-1))
      allocate(npsc(0:numsc-1))

      numsc=0
      do i=0,pnumsc-1
        if(nppsc(i).gt.0) then
          dpsc(i)=numsc
          do j=0,NDTREE-1
            npsc(numsc)=0
            numsc=numsc+1
          enddo
        endif
      enddo

      deallocate(nppsc)

      if(maxlevel.eq.0) then
        if(nscnf.gt.0.and.numsc*NDTREE.lt.MAXNODE) then
! goto lower level
          goto 77
        endif
! *** get the max level among the processors ***
        call MPI_ALLREDUCE(level,maxlevel,1,MPI_INTEGER &
           ,MPI_MAX,MPI_COMM_WORLD,ierr)
        maxlevel = maxlevel+level_mergin
      endif  

      if(level.le.maxlevel) then
        goto 77
      endif

      deallocate(c_nfp)
      deallocate(snfp)
      deallocate(dpsc)
      deallocate(npsc)
      deallocate(scp)

      allocate(phkey1(1:np))
      allocate(ord1(1:np))

      do i=0,np-1
        phkey1(i+1)=phkey(i)
      enddo
      call indexxl(np,phkey1,ord1)

      deallocate(phkey1)

      allocate(list(0:np-1))
      allocate(ord(0:np-1))

      do i=0,np-1
        list(i)=ord1(i+1)-1
      enddo
      do i=0,np-1
        ord(list(i))=i
      enddo

!      write(fileo,'(a4,i3.3)') 'ord1',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        write(60,'(3(1pE13.5),4I20)') x_p(i),y_p(i),z_p(i) &
!         ,id_p(i),phkey(i),ord(i),list(i)
!      enddo
!      close(60)

      deallocate(ord1)
      allocate(phkeyr(0:np))
      allocate(tivr(0:np))

      do i=0,np-1
        phkeyr(i)=phkey(i)
        tivr(i)=list(i)
      enddo
      npj=np
! *** define rank for sending and rank for recieving ***

      allocate(srank(0:0))
      allocate(rrank(0:0))

      srank(0) = myrank+1
      if(srank(0).ge.nprocs) then
        srank(0) = srank(0)-nprocs
      endif
      rrank(0) = myrank-1
      if(rrank(0).lt.0) then
        rrank(0)=rrank(0)+nprocs
      endif
! ** start ordering @ other PE ** 

      allocate(istatus(MPI_STATUS_SIZE))
      allocate(istatusr(MPI_STATUS_SIZE))

      do j = 0,nprocs-1
        if(j.ne.0) then
          crank = myrank-j
          if(crank.lt.0) then
            crank = crank+nprocs
          endif
! *** new ordering using the sorted list ***
          nin=0
          do k=0,np-1
            do i=nin,npj-1
              if(phkeyr(tivr(i)).gt.phkey(list(k))) then
                goto 98
              else if(phkeyr(tivr(i)).eq.phkey(list(k)) &
               .and.crank.ge.myrank) then 
                goto 98
              endif
            enddo
   98       ord(list(k))=ord(list(k))+i
            nin=i
          enddo
        endif
! ** Message Passing **
        isend = npj

        allocate(tivs(0:npj))                
        allocate(phkeys(0:npj))                

        do i=0,isend-1
          tivs(i) = tivr(i)
          phkeys(i) = phkeyr(i)
        enddo
  
        deallocate(tivr)
        deallocate(phkeyr)

!        write(6,*) j,' myrank,isend,npj=',myrank,isend,npj

! get npj
        call MPI_ISEND(isend,1,MPI_INTEGER,srank(0),1,MPI_COMM_WORLD &
         ,ireqs(0),ierr)
        call MPI_IRECV(npj,1,MPI_INTEGER,rrank(0),1,MPI_COMM_WORLD &
         ,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)

!        write(6,*) ' myrank,isend,npj=',myrank,isend,npj

        if(j.ne.nprocs-1) then

          allocate(tivr(0:npj))
          allocate(phkeyr(0:npj))

! ** for order **
          call MPI_ISEND(tivs,isend,MPI_INTEGER,srank(0),2 &
           ,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tivr,npj,MPI_INTEGER,rrank(0),2 &
           ,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
! ** for phkey **
          call MPI_ISEND(phkeys,isend,MPI_INTEGER8,srank(0),3 &
           ,MPI_COMM_WORLD,ireqs(2),ierr)
          call MPI_IRECV(phkeyr,npj,MPI_INTEGER8,rrank(0),3 &
           ,MPI_COMM_WORLD,ireqr(2),ierr)
          call MPI_WAIT(ireqs(2),istatus,ierr)
          call MPI_WAIT(ireqr(2),istatusr,ierr)

!          write(fileo,'(a6,i3.3,i3.3)') 'phkeyr',myrank,rrank(0)
!          open(60,file=fileo,status='unknown')
!          do i=0,npj-1
!            write(60,'(2I20)') phkeyr(i),tivr(i)
!           enddo
!           close(60)

        endif

        deallocate(tivs)
        deallocate(phkeys)

      enddo

!      write(fileo,'(a3,i3.3)') 'ord',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        write(60,'(3(1pE13.5),3I20)') x_p(i),y_p(i),z_p(i) &
!         ,id_p(i),phkey(i),ord(i)
!      enddo
!      close(60)
!      stop

      deallocate(phkey)

      deallocate(srank)
      deallocate(rrank)

! ** send and recieve data ***
! ** define srank and rrank **

      allocate(srank(0:nprocs))
      allocate(rrank(0:nprocs))

      do j=0,nprocs-1
        srank(j) = myrank+j
        if(srank(j).ge.nprocs) then
          srank(j) = srank(j)-nprocs
        endif
        rrank(j) = myrank-j
        if(rrank(j).lt.0) then
          rrank(j) = rrank(j)+nprocs
        endif
      enddo
! ** initialization **

      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))

      do i=0,nprocs-1
        call para_range(0,npt-1,nprocs,i,jsta,jend)
        jjlen(i)=jend-jsta+1
        idisp(i)=jsta
!        write(*,*) ' @ ',myrank,' jsta,jend=',jsta,jend        
      enddo
      npj = np

      allocate(jstap(0:nprocs))
      allocate(jendp(0:nprocs))
      allocate(npsp(0:nprocs))
      allocate(nprp(0:nprocs))

! ** get start and end id, and number of particles ***
      do j = 0,nprocs-1
! ** define array and number sended to srank **
        jstap(j) = 0
        do i = 0,npj-1
          if(ord(list(i)).lt.idisp(srank(j))) then
            jstap(j) = i+1
          else if(ord(list(i)).gt.idisp(srank(j))+jjlen(srank(j))-1) then
            jendp(j) = i-1
            goto 999
          endif
        enddo
        jendp(j) = npj-1
! *** particle number to send ***
  999   npsp(j) = jendp(j)-jstap(j)+1
        nprp(j) = 0
      enddo
! *** start data transfer ***
! get number of particles -> nprp and np
      np=0
      do j = 0,nprocs-1
        call MPI_ISEND(npsp(j),1,MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(nprp(j),1,MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        np=np+nprp(j)
!        write(6,*) ' myrank,j,npsp,nprp=',myrank,j,npsp(j),nprp(j)
      enddo

! *** ord ***

! size of the original particles
      allocate(tivs(0:npj))

      do i=0,npj-1
        tivs(i)=ord(list(i))
      enddo

      deallocate(ord)
! size with the new set of particles
      allocate(ord(0:np))

      np = 0
      do i=jstap(0),jendp(0)
        ord(np)=tivs(i)
        np=np+1
      enddo
      do j = 1,nprocs-1
        call MPI_ISEND(tivs(jstap(j)),npsp(j),MPI_INTEGER,srank(j) &
         ,1,MPI_COMM_WORLD,ireqs(0),ierr)
        call MPI_IRECV(ord(np),nprp(j),MPI_INTEGER,rrank(j) &
         ,1,MPI_COMM_WORLD,ireqr(0),ierr)
        call MPI_WAIT(ireqs(0),istatus,ierr)
        call MPI_WAIT(ireqr(0),istatusr,ierr)
        np=np+nprp(j)
      enddo

      deallocate(tivs)

! *** sort ***
      do i=0,np-1
        ord(i)=ord(i)-idisp(srank(0))
      enddo

! *** for integer values ***
      nval=2

      allocate(tivr(0:np*nval))

      np=0
      do ip = 0,nprocs-1
        npj=0

        allocate(tivs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tivs(npj)=flagc_p(list(i))
          tivs(npj+npsp(ip))=id_p(list(i))
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tivs,npsp(ip)*nval,MPI_INTEGER &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tivr(np*nval),nprp(ip)*nval &
           ,MPI_INTEGER,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tivr(i+nprp(ip)*j)=tivs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tivs)

        np=np+nprp(ip)
      enddo

! reallocate integer values
      call reallocate_baryon_int(np)

      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          flagc_p(ord(npj))=tivr(i)
          id_p(ord(npj))=tivr(i+nprp(ip))
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tivr)
  
      nval=11

      allocate(tdvr(0:np*nval))

      np=0
      do ip = 0,nprocs-1
        npj=0

        allocate(tdvs(0:npsp(ip)*nval))

        do i=jstap(ip),jendp(ip)
          tdvs(npj)=x_p(list(i))
          tdvs(npj+npsp(ip))=y_p(list(i))
          tdvs(npj+npsp(ip)*2)=z_p(list(i))
          tdvs(npj+npsp(ip)*3)=vx_p(list(i))
          tdvs(npj+npsp(ip)*4)=vy_p(list(i))
          tdvs(npj+npsp(ip)*5)=vz_p(list(i))
          tdvs(npj+npsp(ip)*6)=m_p(list(i))
          tdvs(npj+npsp(ip)*7)=u_p(list(i))
          tdvs(npj+npsp(ip)*8)=h_p(list(i))
          tdvs(npj+npsp(ip)*9)=rho_p(list(i))
          tdvs(npj+npsp(ip)*10)=myu_p(list(i))
          npj=npj+1
        enddo
        if(ip.gt.0) then
          call MPI_ISEND(tdvs,npsp(ip)*nval,MPI_DOUBLE_PRECISION &
           ,srank(ip),1,MPI_COMM_WORLD,ireqs(1),ierr)
          call MPI_IRECV(tdvr(np*nval),nprp(ip)*nval &
           ,MPI_DOUBLE_PRECISION,rrank(ip),1,MPI_COMM_WORLD,ireqr(1),ierr)
          call MPI_WAIT(ireqs(1),istatus,ierr)
          call MPI_WAIT(ireqr(1),istatusr,ierr)
        else
          nprp(ip)=npsp(ip)
          do j=0,nval-1         
            do i=0,nprp(ip)-1
              tdvr(i+nprp(ip)*j)=tdvs(i+npsp(ip)*j)
            enddo
          enddo
        endif

        deallocate(tdvs)

        np=np+nprp(ip)
      enddo

! reallocate memory space
      call reallocate_baryon_d1(np)

      npj=0
      np=0
      do ip=0,nprocs-1
        do i=np*nval,np*nval+nprp(ip)-1
          x_p(ord(npj))=tdvr(i)
          y_p(ord(npj))=tdvr(i+nprp(ip))
          z_p(ord(npj))=tdvr(i+nprp(ip)*2)
          vx_p(ord(npj))=tdvr(i+nprp(ip)*3)
          vy_p(ord(npj))=tdvr(i+nprp(ip)*4)
          vz_p(ord(npj))=tdvr(i+nprp(ip)*5)
          m_p(ord(npj))=tdvr(i+nprp(ip)*6)
          u_p(ord(npj))=tdvr(i+nprp(ip)*7)
          h_p(ord(npj))=tdvr(i+nprp(ip)*8)
          rho_p(ord(npj))=tdvr(i+nprp(ip)*9)
          myu_p(ord(npj))=tdvr(i+nprp(ip)*10)
          npj=npj+1
        enddo
        np=np+nprp(ip)
      enddo

      deallocate(tdvr)

      if(np.ne.npj) then
        write(6,*) ' Error in ddecb np,npj=',np,npj
        stop
      endif

      deallocate(srank)
      deallocate(rrank)
      deallocate(jstap)
      deallocate(jendp)
      deallocate(idisp)
      deallocate(jjlen)
      deallocate(npsp)
      deallocate(nprp)

      deallocate(istatus)
      deallocate(istatusr)

! ** after care **

      deallocate(ord)
      deallocate(list)
      allocate(list(0:np))

! reset list_ap, but not take into account star and gas difference
      do i=0,np-1
        list_ap(i)=i
      enddo
     

!      write(fileo,'(a4,i3.3)') 'decp',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        pn=list_ap(i)
!        pn=i
!        write(60,'(5(1pE13.5),I10)') x_p(pn),y_p(pn),z_p(pn) &
!        ,m_p(pn),h_p(pn),id_p(pn)
!      enddo
!      close(60)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      write(6,*) 'after decb ng,nag,myrank=',ng,nag,ns,nas,nagravg,myrank

!      stop


end subroutine ddecb

