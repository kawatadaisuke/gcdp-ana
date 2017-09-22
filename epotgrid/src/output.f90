! ****************************************************
!   output.f for epotgrid pver3
!  1 Aug. 2016    written by D.KAWATA
! ****************************************************

subroutine output(step,ngp,flaggp,ng1,ng2,ng3,ng3i,ran1,ran2,ran3)
      use gcdp_const
      use gcdp_system
      use grid_particle

      implicit none
      include 'mpif.h'

      integer,intent(in) :: step,flaggp,ng1,ng2,ng3,ng3i,ngp
      double precision,intent(in) :: ran1(0:1),ran2(0:1),ran3(0:1)

      integer i,j,nval,ic,ip,npj,ngpt
      integer ir,ith,iz
      character fileo*60
! for work
      integer ierr,ntht
      integer istatus(MPI_STATUS_SIZE)
      double precision rp,frp,lnri,lnro,dlnr,phip
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:)
      double precision,allocatable :: tdvs(:),tdvr(:) 
      double precision,allocatable :: vcirc(:),rpr(:)
 

! analyse circular velocity
      if(flaggp.eq.0) then
        allocate(vcirc(0:ng1-1))
        allocate(rpr(0:ng1-1))

        lnri=dlog(ran1(0))
        lnro=dlog(ran1(1))
        if(ng1.gt.2) then
          dlnr=(lnro-lnri)/dble(ng1-2)
        else
          dlnr=0.0d0
        endif

        vcirc(0)=0.0d0
        rpr(0)=0.0d0 
        if(myrank.eq.0) then
          open(60,file='test.dat')
        endif
        do ir=1,ng1-1
! ir=1 grid point at r=rran(0)
          rpr(ir)=dexp(lnri+dlnr*dble(ir-1))
          vcirc(ir)=0.0d0
          do ith=0,ng3i-1
            ip=pngrid(ir,iz0,ith)
            rp=dsqrt(xp(ip)**2+yp(ip)**2)
            frp=(dvxp(ip)*xp(ip)+dvyp(ip)*yp(ip))/rp
            vcirc(ir)=vcirc(ir)+dsqrt(-frp*rp)
! for test
            if(myrank.eq.0) then
              write(60,'(5(1pE13.5))') rp,dsqrt(-frp*rp)*VUKMS,xp(ip),yp(ip),zp(ip)
            endif
          enddo
        enddo
        if(myrank.eq.0) then
          close(60)
        endif
! sum up 
        allocate(tdvr(0:ng1-1))
        do i=0,ng1-1
          tdvr(i)=0.0d0
        enddo
        call MPI_ALLREDUCE(vcirc,tdvr,ng1,MPI_DOUBLE_PRECISION &
          ,MPI_SUM,MPI_COMM_WORLD,ierr)
        do i=0,ng1-1
          vcirc(i)=tdvr(i)/dble(ng3)
        enddo
      
        deallocate(tdvr)
! output
        write(fileo,'(a12,i6.6,a4)') 'output/vcirc',step,'.bin'
        open(61,file=fileo,status='unknown',form='unformatted')
        write(61) ng1
        write(61) lnri,lnro,dlnr
        write(61) rpr
        write(61) vcirc
        close(61)
        write(fileo,'(a12,i6.6,a4)') 'output/vcirc',step,'.asc'
        open(60,file=fileo,status='unknown')
        write(60,'(a5,I6)') '# nr=',ng1
        do i=0,ng1-1
          write(60,'(2(1pE13.5))') rpr(i),vcirc(i)*VUKMS
        enddo
        close(60)

        deallocate(rpr)
        deallocate(vcirc)

      endif

! for grid output
      nval=7
      ngpt=ng1*ng2*ng3

      allocate(tdvs(0:nval*ngp-1))
      allocate(npjr(0:nprocs-1))
      allocate(ireqs(0:nprocs-1))
      allocate(ireqr(0:nprocs-1))

      if(myrank.eq.0) then
        if(flaggp.eq.0) then
          write(fileo,'(a15,i6.6)') 'output/epotgrid',step
        else
          write(fileo,'(a17,i6.6)') 'output/epotsqgrid',step
        endif
        write(6,*) 'output file name=',fileo
        open(61,file=fileo,status='unknown',form='unformatted')
! nth or ny =ng3
        write(61) ng1,ng3,ng2
        write(61) ran1(0),ran1(1)
        if(flaggp.eq.0) then
          write(61) ran2(0),ran2(1)
        else
          write(61) ran3(0),ran3(1)
          write(61) ran2(0),ran2(1)
        endif
! for test
        if(flaggp.eq.0) then
          write(fileo,'(a12,i6.6)') 'output/aepot',step
        else
          write(fileo,'(a14,i6.6)') 'output/aepotsq',step
        endif
        write(6,*) 'output file name=',fileo
        open(60,file=fileo,status='unknown')
        if(flaggp.eq.0) then
          write(60,*) '# nr,nth,nz=',ng1,ng3,ng2
          write(60,*) '# r range=',ran1(0),ran1(1)
          write(60,*) '# z range=',ran2(0),ran2(1)
        else
          write(60,*) '# nx,ny,nz=',ng1,ng3,ng2
          write(60,*) '# x range=',ran1(0),ran1(1)
          write(60,*) '# y range=',ran3(0),ran3(1)
          write(60,*) '# z range=',ran2(0),ran2(1)
        endif
      endif

! set the output data
      do i=0,ngp-1
        tdvs(i)=xp(i)
        tdvs(i+ngp)=yp(i)
        tdvs(i+ngp*2)=zp(i)
        tdvs(i+ngp*3)=epotp(i)
        tdvs(i+ngp*4)=dvxp(i)
        tdvs(i+ngp*5)=dvyp(i)
        tdvs(i+ngp*6)=dvzp(i)
      enddo

! *** set number of receiving data 
      if(myrank.eq.0) then
        do ip=0,nprocs-1
          if(ip.eq.myrank) then
            npjr(ip)=ngp
          else
            npjr(ip)=0
            call MPI_IRECV(npjr(ip),1,MPI_INTEGER,ip,ip &
              ,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
        enddo
      else
        ip=myrank
        call MPI_ISEND(ngp,1,MPI_INTEGER,0,myrank,MPI_COMM_WORLD,ireqs(ip),ierr)
        call MPI_WAIT(ireqs(ip),istatus,ierr)
      endif    

      if(myrank.eq.0) then
        do ip=0,nprocs-1

          allocate(tdvr(0:npjr(ip)*nval-1))

          if(ip.eq.myrank) then
            do i=0,npjr(myrank)*nval-1
              tdvr(i)=tdvs(i)
            enddo
          else
            call MPI_IRECV(tdvr,npjr(ip)*nval,MPI_DOUBLE_PRECISION &
             ,ip,ip,MPI_COMM_WORLD,ireqr(ip),ierr)
            call MPI_WAIT(ireqr(ip),istatus,ierr)
          endif
! *** write the data
          npj=npjr(ip)
          do i=0,npj-1
            write(61) tdvr(i),tdvr(i+npj),tdvr(i+npj*2),tdvr(i+npj*3) &
              ,tdvr(i+npj*4),tdvr(i+npj*5),tdvr(i+npj*6)
            rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2)
            frp=(tdvr(i+npj*4)*tdvr(i)+tdvr(i+npj*5)*tdvr(i+npj))/rp
            if(rp.gt.0.0d0) then
! angle from -x,y=0, anti-clockwise, y<0, negative phi
              phip=acos(tdvr(i)/rp)
              phip=M_PI-phip
              if(tdvr(i+npj).lt.0.0) then
                phip=-phip
              endif
            else 
              phip=0.0d0
            endif
            write(60,'(10(1pE13.5))') tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
             ,tdvr(i+npj*3),rp*LUKPC,dsqrt(-frp*rp)*VUKMS &
             ,tdvr(i+npj*4),tdvr(i+npj*5),tdvr(i+npj*6),phip
          enddo

          deallocate(tdvr)

        enddo
        close(61)
        close(60)
      else
        ip=myrank
        call MPI_ISEND(tdvs,ngp*nval,MPI_DOUBLE_PRECISION,0,ip &
          ,MPI_COMM_WORLD,ireqs(ip),ierr)
        call MPI_WAIT(ireqs(ip),istatus,ierr)
      endif    

      deallocate(tdvs)
      deallocate(npjr)
      deallocate(ireqs)
      deallocate(ireqr)


end subroutine
