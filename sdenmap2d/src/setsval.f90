! ****************************************************
!    setsval.f90 for sdenmap2d
!  28  Apr. 2015    written by D.KAWATA
! ****************************************************

subroutine setsval(np,nx,ny,xrange,yrange)
      use gcdp_const
      use gcdp_system
      use particle
      use gcdp_kernel
      use mesh

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np,nx,ny
      double precision,intent(in) :: xrange(0:1),yrange(0:1)

! parallerise along x-axis
      integer ixs,ixe
      double precision dx,dy,xmesh,ymesh
! for work
      integer i,j,pn,ierr,ip
      integer is
      double precision xij,yij,rij,hsi,ihsi,s,wij,sdencgs
! mesh data point
      integer ixm0,ixm1,iym0,iym1
      double precision xsm,xem,ysm,yem
      double precision xpm,ypm,zpm
! for work
      integer ntm,nval,nc,nct
      double precision,allocatable :: ty(:),tz(:)
      double precision,allocatable :: tdvs(:),tdvr(:)
      character filen*60

! generate mesh
      allocate(x_m(0:nx-1,0:ny-1))
      allocate(y_m(0:nx-1,0:ny-1))
      allocate(denxy_m(0:nx-1,0:ny-1))
      allocate(denxz_m(0:nx-1,0:ny-1))
      allocate(vrot_m(0:nx-1,0:ny-1))
      allocate(vrad_m(0:nx-1,0:ny-1))
      allocate(vz_m(0:nx-1,0:ny-1))
      allocate(met_m(0:nx-1,0:ny-1))

! assign x,y position
      dx=(xrange(1)-xrange(0))/dble(nx)
      dy=(yrange(1)-yrange(0))/dble(ny)

      do j=0,ny-1
        ypm=yrange(0)+0.5d0*dy+dble(j)*dy
        do i=0,nx-1
          xpm=xrange(0)+0.5d0*dx+dble(i)*dx
          x_m(i,j)=xpm
          y_m(i,j)=ypm
! initialisation
          denxy_m(i,j)=0.0d0
          denxz_m(i,j)=0.0d0
          vrot_m(i,j)=0.0d0
          vrad_m(i,j)=0.0d0
          vz_m(i,j)=0.0d0
          met_m(i,j)=0.0d0
        enddo
      enddo

! at the start and end point of mesh for x and y
      xsm=x_m(0,0)
      xem=x_m(nx-1,0)
      ysm=y_m(0,0)
      yem=y_m(0,ny-1)
! for test
!      write(filen,'(a5,i4.4)') 'contp',myrank
!      open(60,file=filen,status='unknown')

! ip=0: x-y, ip=1: x-z
      do ip=0,1
! swap y and z
        if(ip.eq.1) then

          allocate(ty(0:np-1))
          allocate(tz(0:np-1))

          do i=0,np-1
            ty(i)=yp(i)           
            tz(i)=zp(i)
            yp(i)=tz(i)
            zp(i)=ty(i)
          enddo
        endif

        nc=0
        do pn=0,np-1
! exclude particles obviously not contributing
          if(yp(pn).gt.ysm-hp(pn).and.yp(pn).lt.yem+hp(pn)) then
            if(xp(pn).gt.xsm-hp(pn).and.xp(pn).lt.xem+hp(pn)) then
              nc=nc+1
! output contributing particles
!            write(60,'(4(1pE13.5))') xp(pn),yp(pn),zp(pn),massp(pn)
! x
              ixm0=(xp(pn)-hp(pn)-xsm)/dx
              if(ixm0.lt.0) then
                ixm0=0
              else if(ixm0.ge.nx) then
                ixm0=nx-1
              endif
              ixm1=(xp(pn)+hp(pn)-xsm)/dx
              if(ixm1.lt.0) then
                ixm1=0
              else if(ixm1.ge.nx) then
                ixm1=nx-1
              endif
! y
              iym0=(yp(pn)-hp(pn)-ysm)/dy
              if(iym0.lt.0) then
                iym0=0
              else if(iym1.ge.ny) then
                iym0=ny-1
              endif
              iym1=(yp(pn)+hp(pn)-y_m(0,0))/dy
              if(iym1.lt.0) then
                iym1=0
              else if(iym1.ge.ny) then
                iym1=ny-1
              endif
              do j=iym0,iym1
                do i=ixm0,ixm1
! *** smooth the values ***
                  xij = xp(pn)-x_m(i,j)
                  yij = yp(pn)-y_m(i,j)
                  rij=dsqrt(xij**2+yij**2)
                  hsi=hp(pn)
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
                    if(ip.eq.0) then
                      denxy_m(i,j)=denxy_m(i,j)+massp(pn)*wij
! velocity field only for face-on
                      vrot_m(i,j)=vrot_m(i,j)+massp(pn)*vrotp(pn)*wij
                      vrad_m(i,j)=vrad_m(i,j)+massp(pn)*vradp(pn)*wij
                      vz_m(i,j)=vz_m(i,j)+massp(pn)*vzp(pn)*wij
                      met_m(i,j)=met_m(i,j)+massp(pn)*metp(pn)*wij
                    else
                      denxz_m(i,j)=denxz_m(i,j)+massp(pn)*wij
                    endif
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
! bring back y and z
        if(ip.eq.1) then
          do i=0,np-1
            yp(i)=ty(i)
            zp(i)=tz(i)
          enddo

          deallocate(ty)
          deallocate(tz)

        endif

      enddo
!      close(60)

! sum up
      if(nprocs.gt.1) then
        ntm=nx*ny
        allocate(tdvs(0:ntm-1))
        allocate(tdvr(0:ntm-1))
        nval=6
        do ip=0,nval-1
          nc=0
          if(ip.eq.0) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=denxy_m(i,j)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.1) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=denxz_m(i,j)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.2) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=vrot_m(i,j)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.3) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=vrad_m(i,j)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.4) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=vz_m(i,j)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.5) then
            do j=0,ny-1
              do i=0,nx-1
                tdvs(nc)=met_m(i,j)
                nc=nc+1
              enddo
            enddo
          endif
          do i=0,ntm-1
            tdvr(i)=0.0d0
          enddo
          call MPI_ALLREDUCE(tdvs,tdvr,ntm,MPI_DOUBLE_PRECISION &
           ,MPI_SUM,MPI_COMM_WORLD,ierr)
          nc=0
          if(ip.eq.0) then
            do j=0,ny-1
              do i=0,nx-1
                denxy_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.1) then
            do j=0,ny-1
              do i=0,nx-1
                denxz_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.2) then
            do j=0,ny-1
              do i=0,nx-1
                vrot_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.3) then
            do j=0,ny-1
              do i=0,nx-1
                vrad_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.4) then
            do j=0,ny-1
              do i=0,nx-1
                vz_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          else if(ip.eq.5) then
            do j=0,ny-1
              do i=0,nx-1
                met_m(i,j)=tdvr(nc)
                nc=nc+1
              enddo
            enddo
          endif
        enddo
        deallocate(tdvs)
        deallocate(tdvr)
      endif       

! get the mean value
      do j=0,ny-1
        do i=0,nx-1
          if(denxy_m(i,j).gt.0.0d0) then
            vrot_m(i,j)=vrot_m(i,j)/denxy_m(i,j)
            vrad_m(i,j)=vrad_m(i,j)/denxy_m(i,j)
            vz_m(i,j)=vz_m(i,j)/denxy_m(i,j)
            met_m(i,j)=met_m(i,j)/denxy_m(i,j)
          endif
        enddo
      enddo

end subroutine
      
