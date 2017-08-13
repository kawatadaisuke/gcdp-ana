! ****************************************************
!    setsval.f90 for smadrzmet
!  12  Feb. 2016    written by D.KAWATA
! ****************************************************

subroutine setsval(np,xp,yp,hm,nx,ny,xl,xh,yl,yh,fval)
      use gcdp_const
      use gcdp_system
      use gcdp_kernel

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np,nx,ny
      double precision,intent(in) :: hm,xl,xh,yl,yh
      double precision,intent(in) :: xp(0:np-1),yp(0:np-1)
      double precision,intent(inout) :: fval(0:nx-1,0:ny-1)

      integer ierr
! parallerise along x-axis
      integer ixs,ixe
      double precision dx,dy,xmesh,ymesh
! for work
      integer i,j,pn
      integer is
      double precision xij,yij,rij,hsi,ihsi,s,wij,sdencgs
! mesh data point
      integer ixm0,ixm1,iym0,iym1
      double precision xsm,xem,ysm,yem
      double precision xpm,ypm,zpm
! for MPI work
      integer ntm,nc,nct
      double precision,allocatable :: tdvs(:),tdvr(:)

!     write(6,*) ' xl,xh,yl,yh=',xl,xh,yl,yh

! assign x,y position
      dx=(xh-xl)/dble(nx)
      dy=(yh-yl)/dble(ny)

! at the start and end point of mesh for x and y
      xsm=xl+0.5d0*dx
      xem=xh-0.5d0*dx
      ysm=yl+0.5d0*dy
      yem=yh-0.5d0*dy

      nc=0
      if(myrank.eq.0) then
        write(6,*) ' np in setsval=',np
      endif
!        write(6,*) myrank,hm,xl,xh,yl,yh
      do pn=0,np-1
! exclude particles obviously not contributing
        if(yp(pn).gt.ysm-hm.and.yp(pn).lt.yem+hm) then
          if(xp(pn).gt.xsm-hm.and.xp(pn).lt.xem+hm) then
            nc=nc+1
! output contributing particles
!            write(60,'(2(1pE13.5))') xp(pn),yp(pn)
! x
            ixm0=(xp(pn)-hm-xsm)/dx
            if(ixm0.lt.0) then
              ixm0=0
            else if(ixm0.gt.nx-1) then
              ixm0=nx-1
            endif
            ixm1=(xp(pn)+hm-xsm)/dx
            if(ixm1.lt.0) then
              ixm1=0
            else if(ixm1.gt.nx-1) then
              ixm1=nx-1
            endif
! y
            iym0=(yp(pn)-hm-ysm)/dy
            if(iym0.lt.0) then
              iym0=0
            else if(iym1.gt.ny-1) then
              iym0=ny-1
            endif
            iym1=(yp(pn)+hm-ysm)/dy
            if(iym1.lt.0) then
              iym1=0
            else if(iym1.gt.ny-1) then
              iym1=ny-1
            endif
            do j=iym0,iym1
              ypm=yl+0.5d0*dy+dble(j)*dy
              do i=ixm0,ixm1
                xpm=xl+0.5d0*dx+dble(i)*dx
! *** smooth the values ***
                xij = xp(pn)-xpm
                yij = yp(pn)-ypm
                rij=dsqrt(xij**2+yij**2)
                hsi=hm
                if(rij.lt.hsi) then
                  ihsi=1.0d0/hsi
                  s=rij*ihsi
                  is=int(s*dnktab)
                  if(is.lt.0) then
                    is=0
                  else if(is.ge.NKTAB) then
                    is=NKTAB-1
                  endif
                  wij=w_tb(is)+(w_tb(is+1)-w_tb(is))*(s-s_tb(is))*dnktab
! 2D
                  wij=wij*(ihsi**2)
                  fval(i,j)=fval(i,j)+wij
                endif
              enddo
            enddo
          endif
        endif
      enddo
      nct=0
      call MPI_ALLREDUCE(nc,nct,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0) then
        write(6,*) ' # of contributing particles=',nct
      endif
!      close(60)

! sum up
      if(nprocs.gt.1) then
        ntm=nx*ny
        allocate(tdvs(0:ntm-1))
        allocate(tdvr(0:ntm-1))

        nc=0
        do j=0,ny-1
          do i=0,nx-1
            tdvs(nc)=fval(i,j)
            tdvr(nc)=0.0d0
            nc=nc+1
          enddo
        enddo
        call MPI_ALLREDUCE(tdvs,tdvr,ntm,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,MPI_COMM_WORLD,ierr)
        nc=0
        do j=0,ny-1
          do i=0,nx-1
            fval(i,j)=tdvr(nc)
            nc=nc+1
          enddo
        enddo

        deallocate(tdvr)
        deallocate(tdvs)

      endif

end subroutine
      
