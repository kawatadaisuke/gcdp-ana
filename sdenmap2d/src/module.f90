! ****************************************
!   module.f90 for sden2d
!  18 Jan. 2016  written by D.KAWATA
! ****************************************

module gcdp_const
      implicit none
      
      double precision,parameter :: M_PI=3.141592653589793238d0

! *** for Tree ***
! number of dauthter tree
      integer,parameter :: NDTREE=8
! maximum number of node for ddec
      integer,parameter :: MAXNODE=80000000
! maximum number of node to be sent
      integer,parameter :: MAXNODESEND=70
! Margin for Size of root tree
      double precision,parameter :: MGROOT=1.001d0
! Max number of particles
      integer,parameter :: MNB=10000000

! *** for at Table ***
      integer,parameter :: NATTABLE=50000
! *** for kernel tables ***
      integer,parameter :: NKTAB=100000

! *** For Arithmetic Value ***
      double precision,parameter :: THIRD=1.0d0/3.0d0
      double precision,parameter :: V34=3.0d0/4.0d0
      double precision,parameter :: INF=1.0d6
      double precision,parameter :: MININF=1.0d-20
      double precision,parameter :: MININFESN=1.0d-20

! *** Normalization Unit ***
! * mass unit ( solar mass ) *      
      double precision,parameter :: MUSM=1.0d12
! * (solarmass)/(gram) *      
!      parameter (MUG=1.99e45)
! * length unit / pc ( 100kpc/pc ) *      
      double precision,parameter :: LUPC=1.0e5
! * length unit(100kpc)/kpc *      
      double precision,parameter :: LUKPC=1.0e2
! * (100 kpc)/(cm) *      
      double precision,parameter :: LUCM=3.086e23
! * Density *      
      double precision,parameter :: DU=6.77e-26
! * Pressure *      
      double precision,parameter :: PU=2.91e-11
! * tempature ( K ) *      
      double precision,parameter :: TUK=1.0e4
! * time ( yr ) *      
      double precision,parameter :: TMUYR=4.71e8
! * time ( Gyr ) *      
      double precision,parameter :: TMUGYR=0.471d0
! * time (s) *      
      double precision,parameter :: TMUS=1.488d16
! * km/s *      
      double precision,parameter :: VUKMS=207.4d0
! * cooling (energy) rate (erg s^-1 cm^-3) *      
      double precision,parameter :: ERU=1.96d-27
! * k/m unit /( cm^2 s^-2 K^-1) *       
      double precision,parameter :: K_MU=4.3d10

! *** For Physical Value ***
! * Gravitatrional Constant *      
      double precision,parameter :: G=1.0d0
! * Boltzmann constarnt *    
      double precision,parameter :: KCGS=1.381d-16
! * average molecular weight (default) *      
      double precision,parameter :: MYU=0.6d0
! * proton mass *      
      double precision,parameter :: MP=1.67265d-24
! * G(cgs) *      
      double precision,parameter :: GCGS=6.672d-8
      double precision,parameter :: TPRHO=((KCGS/(MYU*MP))/K_MU)
! * Hubble constant 9.78x10^9xh^-1 yr h = 1.0*      
      double precision,parameter :: H0_1=(1.0d0/((100.0d0*1.0d5)/(3.086d24))) &
       /((((3600.0d0*24.0d0)*365.0d0)*4.0d0+(3600.0d0*24.0d0))/4.0d0)
      double precision,parameter :: HUB0=100.0d0/(VUKMS*10.0d0)
! * Solar metallicity *
      double precision,parameter :: XZSOL=0.019d0
! He primordial and Solar metallicity
      double precision,parameter :: XHE0=0.24d0
      double precision,parameter :: XHESOL=0.275d0
! * photon velocity *
      double precision,parameter :: CVEL=(2.997925e5/VUKMS)

! *** note:  h_gcdp = 2 h_monaghan ***
      double precision,parameter :: ETAH=2.7d0
      double precision,parameter :: ERRH=0.01d0
! for setting stellar h_p
      double precision,parameter :: NSTH=0.1d0
      double precision,parameter :: XHSOL=0.706d0
! *** limit for the maximum change of h during the iteration ***
      double precision,parameter :: DHFLIM=4.0d0

! * For Parallel *
! maximum number of cpus
      integer,parameter :: NCPU=256
! *** Peano-Hilbert curve ***
      integer*8,parameter :: NDPH=8

end module gcdp_const

module gcdp_baryon

      implicit none
! * Particle ID *
      integer,allocatable,save :: id_p(:)
! * number neighbour *
      integer,allocatable,save :: nnb_p(:)
! * Mass *      
      double precision,allocatable,save :: m_p(:)
! * Virtual Position Xn+1,Xn *      
      double precision,allocatable,save :: x_p(:),y_p(:),z_p(:)
! * Velosity Vn-1/2,Vn+1/2 *
      double precision,allocatable,save :: vx_p(:),vy_p(:),vz_p(:)
      double precision,allocatable,save :: u_p(:)
! * molecular weight *
      double precision,allocatable,save :: myu_p(:)
! * Smoothing Length hn+1,hn *      
      double precision,allocatable,save :: h_p(:)
! * Density rhon+1,rhon *       
      double precision,allocatable,save :: rho_p(:)
! * Presure Pn+1,Pn *      
      double precision,allocatable,save :: p_p(:)
! * Sound Velocity *      
      double precision,allocatable,save :: cs_p(:)
! * Entropic function *
      double precision,allocatable,save :: as_p(:)
! * div V *
      double precision,allocatable,save :: div_v_p(:)
! *** artificial viscosity and thermal conductivity term ***
      double precision,allocatable,save :: alpv_p(:),alpu_p(:)
! ***** For Individual Time Step *****
! * list of active particles for gas and star *      
      integer,allocatable,save :: list_ap(:)

! * For Heavy Elements unit Solar Mass *
      double precision,allocatable,save :: mzHe_p(:),mzC_p(:),mzN_p(:) &
       ,mzO_p(:),mzNe_p(:),mzMg_p(:) &
       ,mzSi_p(:),mzFe_p(:),mzZ_p(:)

      integer,allocatable,save :: flagfd_p(:),flagc_p(:)

end module gcdp_baryon

module gcdp_dm
      use gcdp_const
      implicit none
! ***** For Dark Matter *****
! * Particle ID *
      integer,allocatable,save :: id_dm(:)
! *** number of neighbour particles ***
      integer,allocatable,save :: nnb_dm(:)
! * mass *      
      double precision,allocatable,save :: m_dm(:)
! * Virtual Position *      
      double precision,allocatable,save :: x_dm(:),y_dm(:),z_dm(:)
! * Velocity V(n+1/2) *
      double precision,allocatable,save :: vx_dm(:),vy_dm(:),vz_dm(:)
! *** for adaptive softening ***
      double precision,allocatable,save :: h_dm(:),rho_dm(:)
! * list of active particles for DM *      
      integer,allocatable,save :: list_adm(:)

end module gcdp_dm

module gcdp_system
      implicit none

! ***** for MPI *****
      integer,save :: nprocs,myrank

end module gcdp_system

module gcdp_kernel
      use gcdp_const
      implicit none
! *** kernel look up table ****
      double precision,save :: dnktab
      double precision,save :: s_tb(0:NKTAB),w_tb(0:NKTAB),dwds_s_tb(0:NKTAB)
end module gcdp_kernel

module gcdp_gtree
      use gcdp_const
      implicit none
! ***** Tree for SPH particles _gtr *****
      integer,save :: num_gtr
! * max node id sent to the other proc 
      integer,save :: nodese_gtr,nodess_gtr
! * number of contained particle *       
      integer,allocatable,save :: np_gtr(:)
! * name of Particle *            
      integer,allocatable,save :: pn_gtr(:)
! * length of side *      
      double precision,allocatable,save :: l_gtr(:),hm_gtr(:)
      double precision,allocatable,save :: cx_gtr(:),cy_gtr(:),cz_gtr(:)
! * first child node *      
      integer,allocatable,save :: daughter_gtr(:)
! * next node *      
      integer,allocatable,save :: next_gtr(:)
! *** proc id ***
      integer,allocatable,save :: proc_gtr(:)
end module gcdp_gtree

module particle
      integer,allocatable,save :: idp(:),listp(:)
      double precision,allocatable,save :: xp(:),yp(:),zp(:),massp(:) &
       ,hp(:),vxp(:),vyp(:),vzp(:),metp(:),vrotp(:),vradp(:)
end module particle

module mesh
      double precision,allocatable,save :: x_m(:,:),y_m(:,:)
      double precision,allocatable,save :: denxy_m(:,:),denxz_m(:,:)
      double precision,allocatable,save :: vrot_m(:,:),vrad_m(:,:),vz_m(:,:)
      double precision,allocatable,save :: met_m(:,:)
end module mesh
