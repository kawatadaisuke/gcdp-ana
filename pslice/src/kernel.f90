! *****************************************************
!    kernel.f for 3dgrid
!  14  Aug. 2014   written by D. Kawata
! ***************************************************** 

subroutine setkernel(flag3d)
      use gcdp_const
      use gcdp_kernel
      use gcdp_system

      implicit none

      integer,intent(in) :: flag3d
      integer i,j
      double precision f,cv,s,ds_tb

      dnktab=dble(NKTAB)
      ds_tb=1.0d0/dnktab

! *** W: SPH kernel ***
      if (flag3d.le.0) then
        cv=8.0d0/(M_PI)
      else
! *** for 2D constant ***
        cv=(40.0d0/(7.0d0*M_PI))         
      endif

      do i=0,NKTAB
        s_tb(i)=ds_tb*dble(i)
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(1.0d0-6.0d0*(s**2)+6.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(2.0d0*((1.0d0-s)**3))
        else 
          f=0.0d0
        endif
        w_tb(i)=cv*f
      enddo

! *** dW(s)/ds/s: SPH kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=-12.0d0+18.0d0*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)/s
        else 
          f=0.0d0
        endif
        dwds_s_tb(i)=cv*f
      enddo

end subroutine
