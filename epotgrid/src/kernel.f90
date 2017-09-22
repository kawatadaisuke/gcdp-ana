! *****************************************************
!    kernel.F95 for GCD+ ver. f03.0
! 04  Jan. 2013   written by D. Kawata
! ***************************************************** 

subroutine setkernel()
      use gcdp_const
      use gcdp_kernel
      use gcdp_system

      implicit none

      integer i,j
      double precision f,cv,s,ds_tb

      dnktab=dble(NKTAB)
      ds_tb=1.0d0/dnktab

! *** W: SPH kernel ***
! *** for constant ***
      cv=8.0d0/M_PI

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
! *** core dw/ds ***
        if(s.le.THIRD) then
          f=-2.0d0
        else if(s.le.0.5d0) then
          f=-12.0d0*s+18.0d0*s*s
        else if(s.le.1.0d0) then
          f=-6.0d0*((1.0d0-s)**2)
        else 
          f=0.0d0
        endif
        dwdsc_tb(i)=cv*f
      enddo
! *** h x dW/dh: SPH kernel ***
! #ifdef SIM1D
!      cv=1.0d0
! #elif defined(SIM2D)
!      cv=2.0d0
! #else
!      cv=3.0d0
! #endif
!      do i=0,NKTAB
!        s=s_tb(i)
!        hdwdh_tb(i)=-cv*w_tb(i)-dwds_s_tb(i)*(s**2)
!      enddo

! *** h*phi: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.lt.0.5d0) then
          f=(16.0d0/3.0d0)*(s**2)-(48.0d0/5.0d0)*(s**4)+(32.0d0/5.0d0)*(s**5) &
            -(14.0d0/5.0d0)
        else
          f=(32.0d0/3.0d0)*(s**2)-16.0d0*(s**3)+(48.0d0/5.0d0)*(s**4) &
            -(32.0d0/15.0d0)*(s**5)-(16.0d0/5.0d0)+(1.0d0/(15.0d0*s))
        endif
        hphi_tb(i)=f
      enddo

! *** dphi/dr/r: phi kernel ***
      do i=0,NKTAB
        s=s_tb(i)
        if(s.le.0.5d0) then
          f=(32.0d0/15.0d0)*(5.0d0-18.0d0*(s**2)+15.0d0*(s**3))
        else if(s.le.1.0d0) then
          f=(4.0d0/15.0d0)*(80.0d0*s-180.0d0*(s**2)+144.0d0*(s**3) &
           -40.0d0*(s**4)-0.25d0/(s**2))/s
        endif
        dphidr_r_tb(i)=f
      enddo

! *** for test
!      if(myrank.eq.0) then
!        open(60,file='kernel.dat',status='unknown')
!        do i=1,NKTAB-1
!          s=s_tb(i)
!          write(60,160) s_tb(i),w_tb(i),hphi_tb(i)
!           ,dwds_s_tb(i),(dwds_s_tb(i)-dwds_s_tb(i-1))/dwds_s_tb(i) &
!           ,dphidr_r_tb(i),(dphidr_r_tb(i)-dphidr_r_tb(i-1))/dphidr_r_tb(i) &
!           ,dphidh_tb(i),(dphidh_tb(i)-dphidh_tb(i-1))/dphidh_tb(i) &
!           ,d2phidr2_tb(i),(d2phidr2_tb(i)-d2phidr2_tb(i-1))/d2phidr2_tb(i) &
!           ,d3phidr3_tb(i),(d3phidr3_tb(i)-d3phidr3_tb(i-1))/d3phidr3_tb(i)
!        enddo
!  160   format(3(1pE13.5))
!        close(60)
!      endif

end subroutine
