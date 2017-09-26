! ****************************************************
!   output.f for sdenmap3d
!  28 Apr. 2015    written by D.KAWATA
! ****************************************************

subroutine output(step,nxgrid,nygrid,xrange,yrange,tu,flagc,flaginc,iang)
      use gcdp_const
      use gcdp_system
      use mesh

      implicit none
      include 'mpif.h'

      integer,intent(in) :: step,nxgrid,nygrid,flagc,flaginc
      double precision,intent(in) :: xrange(0:1),yrange(0:1),tu
      double precision,intent(inout) :: iang

      integer i,j,nval,ic,nmp,nmop,ip,npj
      character fileo*60
      double precision tym,tzmi,angc
! for work
      integer ierr
      integer istatus(MPI_STATUS_SIZE)
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)    

      if(myrank.eq.0) then
        if(flagc.eq.0) then
          write(fileo,'(a13,i6.6)') 'output/gden2d',step
        else if(flagc.eq.1) then
          write(fileo,'(a13,i6.6)') 'output/dden2d',step
        else
          write(fileo,'(a13,i6.6)') 'output/sden2d',step
        endif
        write(6,*) 'output file name=',fileo
        open(61,file=fileo,status='unknown',form='unformatted')
        write(61) nxgrid,nygrid
        write(61) xrange(0),xrange(1)
        write(61) yrange(0),yrange(1)
        write(61) tu*TMUGYR
        write(61) flaginc
        iang=iang*180.0d0/M_PI
        write(61) iang
! for ascii data
        if(flagc.eq.0) then
          write(fileo,'(a15,i6.6,a4)') 'output/gden2dxy',step,'.asc'
        else if(flagc.eq.1) then
          write(fileo,'(a15,i6.6,a4)') 'output/dden2dxy',step,'.asc'
        else
          write(fileo,'(a15,i6.6,a4)') 'output/sden2dxy',step,'.asc'
        endif
        write(6,*) 'output ascii file name=',fileo
        open(60,file=fileo,status='unknown')
        write(60,'(a1,2I10)') '#',nxgrid,nygrid
        write(60,'(a1,2(1pE13.5))') '#',xrange(0),xrange(1)
        write(60,'(a1,2(1pE13.5))') '#',yrange(0),yrange(1)
        write(60,'(a1,1pE13.5)') '#',tu*TMUGYR
        write(60,'(a1,I3,1pE13.5)') '#',flaginc,iang

        do j=0,nygrid-1
          do i=0,nxgrid-1
            write(61) x_m(i,j),y_m(i,j),denxy_m(i,j),denxz_m(i,j) &
             ,vrot_m(i,j),vrad_m(i,j),vz_m(i,j),met_m(i,j)
          enddo
        enddo

! deg -> radian
        if(flaginc.eq.0) then
          angc=1.0d0
        else
          angc=1.0d0/dcos(iang*M_PI/180.0d0)
        endif
        do j=0,nygrid-1
          do i=0,nxgrid-1
! correct inclination only for y coordinate
            tym=y_m(i,j)*angc
            write(60,'(8(1pE13.5))') x_m(i,j),y_m(i,j) &
             ,denxy_m(i,j)*MUSM/(1.0e6) &
             ,vrot_m(i,j),vrad_m(i,j),vz_m(i,j),met_m(i,j) &
             ,dsqrt(x_m(i,j)**2+tym**2)
          enddo
        enddo
        close(60)
        close(61)

      endif

end subroutine
