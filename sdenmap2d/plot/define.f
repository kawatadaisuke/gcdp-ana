c/****************************************
c   define.f for spp ver.1 
c 2   Oct. '98    produced by D.KAWATA
c ****************************************/

      implicit none

      double precision M_PI
      parameter (M_PI=3.141592654)     
c/*** Definition Max Number ***/
c/*** for particles ***/
c/* Max Number of Particle */
      integer MNR
      parameter (MNR=200)
c/* Max resolutsion in 1 D */
      integer MR
      parameter (MR=1100)

c/*** For Arithmetic Value ***/
      double precision THIRD,INF,DMINV
      parameter (THIRD=1.0/3.0)
      parameter (INF=1000000.0)
      parameter (DMINV=1.0e-6)

c /**** Simulattion Unit *****/
      double precision MUSM,TMUYR,LUKPC,VUKMS,DU
c /* mass unit ( solar mass ) */
      parameter(MUSM=1.0e12)
c /* time ( yr ) */
      parameter (TMUYR=4.71e8)
c /* length unit(100kpc)/kpc */
      parameter (LUKPC=1.0e2)
c /* km/s */
      parameter (VUKMS=207.4d0)
c /* Density */      
      parameter (DU=6.77e-26)

c /*** For correction to massive star ***/
      double precision MSNR,MFSN,NFSN,MFSN11,NFSN11
c /* Mass fraction of Massive star */
      parameter (MSNR=1.4)      
      parameter (MFSN=0.122182)
      parameter (NFSN=0.007310)
      parameter (MFSN11=0.249163)      
      parameter (NFSN11=0.013824)
