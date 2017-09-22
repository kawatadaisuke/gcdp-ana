! *****************************************************
!    allocate.F95 for GCD+ ver. f03.0
! 25  Feb. 2013   written by D. Kawata
! ***************************************************** 

subroutine allocate_baryon(np)
      use gcdp_baryon

      implicit none

      integer,intent(in) :: np

      if(allocated(id_p)) then
        deallocate(id_p)
        deallocate(nnb_p)

        deallocate(m_p)
        deallocate(x_p)
        deallocate(y_p)
        deallocate(z_p)
        deallocate(vx_p)
        deallocate(vy_p)
        deallocate(vz_p)
        deallocate(u_p)
        deallocate(myu_p)
        deallocate(h_p)
        deallocate(rho_p)
        deallocate(alpv_p)
        deallocate(alpu_p)
        deallocate(div_v_p)
        deallocate(list_ap)
        deallocate(flagc_p)
        deallocate(flagfd_p)
      endif

      allocate(id_p(0:np-1))
      allocate(nnb_p(0:np-1))

      allocate(m_p(0:np-1))
      allocate(x_p(0:np-1))
      allocate(y_p(0:np-1))
      allocate(z_p(0:np-1))
      allocate(vx_p(0:np-1))
      allocate(vy_p(0:np-1))
      allocate(vz_p(0:np-1))
      allocate(u_p(0:np-1))
      allocate(myu_p(0:np-1))
      allocate(h_p(0:np-1))
      allocate(rho_p(0:np-1))
      allocate(div_v_p(0:np-1))
      allocate(alpv_p(0:np-1))
      allocate(alpu_p(0:np-1))

      allocate(list_ap(0:np-1))
      allocate(flagc_p(0:np-1))
      allocate(flagfd_p(0:np-1))

end subroutine

! reallocate integer baryon variables
subroutine reallocate_baryon_int(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(id_p)
      deallocate(flagc_p)
      deallocate(list_ap)

! allocate
      allocate(id_p(0:np))
      allocate(flagc_p(0:np))
      allocate(list_ap(0:np))

end subroutine reallocate_baryon_int

! reallocate first set of double precision baryon variables
subroutine reallocate_baryon_d1(np)
      use gcdp_baryon

      implicit none
      integer,intent(in) :: np

! deallocate
      deallocate(x_p)
      deallocate(y_p)
      deallocate(z_p)
      deallocate(vx_p)
      deallocate(vy_p)
      deallocate(vz_p)
      deallocate(m_p)
      deallocate(h_p)
      deallocate(u_p)
      deallocate(rho_p)
      deallocate(myu_p)

! allocate
      allocate(x_p(0:np))
      allocate(y_p(0:np))
      allocate(z_p(0:np))
      allocate(vx_p(0:np))
      allocate(vy_p(0:np))
      allocate(vz_p(0:np))
      allocate(m_p(0:np))
      allocate(u_p(0:np))
      allocate(h_p(0:np))
      allocate(rho_p(0:np))
      allocate(myu_p(0:np))

end subroutine


subroutine allocate_dm(ndm)
      use gcdp_dm

      implicit none

      integer,intent(in) :: ndm

      if(allocated(id_dm)) then
        deallocate(id_dm)
        deallocate(x_dm)
        deallocate(y_dm)
        deallocate(z_dm)
        deallocate(vx_dm)
        deallocate(vy_dm)
        deallocate(vz_dm)
        deallocate(m_dm)
        deallocate(rho_dm)
        deallocate(h_dm)
        deallocate(list_adm)
      endif

      allocate(id_dm(0:ndm-1))
      allocate(x_dm(0:ndm-1))
      allocate(y_dm(0:ndm-1))
      allocate(z_dm(0:ndm-1))
      allocate(vx_dm(0:ndm-1))
      allocate(vy_dm(0:ndm-1))
      allocate(vz_dm(0:ndm-1))
      allocate(m_dm(0:ndm-1))
      allocate(rho_dm(0:ndm-1))
      allocate(h_dm(0:ndm-1))
      allocate(list_adm(0:ndm-1))

end subroutine

! reallocate DM integer variable
subroutine reallocate_dm_int(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(id_dm)
      deallocate(list_adm)

! allocate
      allocate(id_dm(0:ndm))
      allocate(list_adm(0:ndm))

end subroutine

! reallocate DM double variables 1
subroutine reallocate_dm_d1(ndm)
      use gcdp_dm

      implicit none
      integer,intent(in) :: ndm

! deallocate
      deallocate(x_dm)
      deallocate(y_dm)
      deallocate(z_dm)
      deallocate(vx_dm)
      deallocate(vy_dm)
      deallocate(vz_dm)
      deallocate(m_dm)
      deallocate(h_dm)
      deallocate(rho_dm)

! allocate
      allocate(x_dm(0:ndm))
      allocate(y_dm(0:ndm))
      allocate(z_dm(0:ndm))
      allocate(vx_dm(0:ndm))
      allocate(vy_dm(0:ndm))
      allocate(vz_dm(0:ndm))
      allocate(m_dm(0:ndm))
      allocate(h_dm(0:ndm))
      allocate(rho_dm(0:ndm))

end subroutine


! allocate general particle data
subroutine allocate_gridparticle(np,nr,nz,nthi)
      use grid_particle

      implicit none

      integer,intent(in) :: np,nr,nz,nthi

      if(allocated(xp)) then
        deallocate(xp)
        deallocate(yp)
        deallocate(zp)
        deallocate(epotp)
        deallocate(dvxp)
        deallocate(dvyp)
        deallocate(dvzp)
        deallocate(hp)
        deallocate(pngrid)
      endif

      allocate(xp(0:np-1))
      allocate(yp(0:np-1))
      allocate(zp(0:np-1))
      allocate(dvxp(0:np-1))
      allocate(dvyp(0:np-1))
      allocate(dvzp(0:np-1))
      allocate(epotp(0:np-1))
      allocate(hp(0:np-1))
      allocate(pngrid(0:nr-1,0:nz-1,0:nthi-1))

end subroutine

! allocate btree
subroutine allocate_btree(ntr)
      use gcdp_btree

      implicit none
      integer,intent(in) :: ntr

      if(allocated(np_tr)) then
        deallocate(np_tr)
        deallocate(pn_tr)
        deallocate(l_tr)
        deallocate(cx_tr)
        deallocate(cy_tr)
        deallocate(cz_tr)
        deallocate(daughter_tr)
        deallocate(next_tr)
        deallocate(cmx_tr)
        deallocate(cmy_tr)
        deallocate(cmz_tr)
        deallocate(mass_tr)
        deallocate(delta_tr)
        deallocate(hm_tr)
        deallocate(mx_tr)
        deallocate(my_tr)
        deallocate(mz_tr)
        deallocate(mxx_tr)
        deallocate(myy_tr)
        deallocate(mzz_tr)
        deallocate(mxy_tr)
        deallocate(myz_tr)
        deallocate(mzx_tr)
        deallocate(proc_tr)
      endif

      allocate(np_tr(0:ntr))
      allocate(pn_tr(0:ntr))
      allocate(l_tr(0:ntr))
      allocate(cx_tr(0:ntr))
      allocate(cy_tr(0:ntr))
      allocate(cz_tr(0:ntr))
      allocate(daughter_tr(0:ntr))
      allocate(next_tr(0:ntr))
      allocate(cmx_tr(0:ntr))
      allocate(cmy_tr(0:ntr))
      allocate(cmz_tr(0:ntr))
      allocate(mass_tr(0:ntr))
      allocate(delta_tr(0:ntr))
      allocate(hm_tr(0:ntr))
      allocate(mx_tr(0:ntr))
      allocate(my_tr(0:ntr))
      allocate(mz_tr(0:ntr))
      allocate(mxx_tr(0:ntr))
      allocate(myy_tr(0:ntr))
      allocate(mzz_tr(0:ntr))
      allocate(mxy_tr(0:ntr))
      allocate(myz_tr(0:ntr))
      allocate(mzx_tr(0:ntr))
      allocate(proc_tr(0:ntr))

end subroutine

subroutine allocate_dmtree(ntr)
      use gcdp_dmtree

      implicit none
      integer,intent(in) :: ntr

      if(allocated(np_dmtr)) then
        deallocate(np_dmtr)
        deallocate(pn_dmtr)
        deallocate(l_dmtr)
        deallocate(cx_dmtr)
        deallocate(cy_dmtr)
        deallocate(cz_dmtr)
        deallocate(daughter_dmtr)
        deallocate(next_dmtr)
        deallocate(cmx_dmtr)
        deallocate(cmy_dmtr)
        deallocate(cmz_dmtr)
        deallocate(mass_dmtr)
        deallocate(delta_dmtr)
        deallocate(hm_dmtr)

        deallocate(mx_dmtr)
        deallocate(my_dmtr)
        deallocate(mz_dmtr)
        deallocate(mxx_dmtr)
        deallocate(myy_dmtr)
        deallocate(mzz_dmtr)
        deallocate(mxy_dmtr)
        deallocate(myz_dmtr)
        deallocate(mzx_dmtr)

        deallocate(proc_dmtr)
      endif

      allocate(np_dmtr(0:ntr))
      allocate(pn_dmtr(0:ntr))
      allocate(l_dmtr(0:ntr))
      allocate(cx_dmtr(0:ntr))
      allocate(cy_dmtr(0:ntr))
      allocate(cz_dmtr(0:ntr))
      allocate(daughter_dmtr(0:ntr))
      allocate(next_dmtr(0:ntr))
      allocate(cmx_dmtr(0:ntr))
      allocate(cmy_dmtr(0:ntr))
      allocate(cmz_dmtr(0:ntr))
      allocate(mass_dmtr(0:ntr))
      allocate(delta_dmtr(0:ntr))
      allocate(hm_dmtr(0:ntr))

      allocate(mx_dmtr(0:ntr))
      allocate(my_dmtr(0:ntr))
      allocate(mz_dmtr(0:ntr))
      allocate(mxx_dmtr(0:ntr))
      allocate(myy_dmtr(0:ntr))
      allocate(mzz_dmtr(0:ntr))
      allocate(mxy_dmtr(0:ntr))
      allocate(myz_dmtr(0:ntr))
      allocate(mzx_dmtr(0:ntr))

      allocate(proc_dmtr(0:ntr))

end subroutine
