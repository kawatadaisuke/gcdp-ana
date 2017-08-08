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
        deallocate(p_p)
        deallocate(cs_p)
        deallocate(as_p)
        deallocate(div_v_p)
        deallocate(alpv_p)
        deallocate(alpu_p)
        deallocate(list_ap)
        deallocate(mzHe_p)
        deallocate(mzC_p)
        deallocate(mzN_p)
        deallocate(mzO_p)
        deallocate(mzNe_p)
        deallocate(mzMg_p)
        deallocate(mzSi_p)
        deallocate(mzFe_p)
        deallocate(mzZ_p)
        deallocate(ts_p)
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
      allocate(p_p(0:np-1))
      allocate(cs_p(0:np-1))
      allocate(as_p(0:np-1))
      allocate(div_v_p(0:np-1))
      allocate(alpv_p(0:np-1))
      allocate(alpu_p(0:np-1))
      allocate(list_ap(0:np-1))
      allocate(mzHe_p(0:np-1))
      allocate(mzC_p(0:np-1))
      allocate(mzN_p(0:np-1))
      allocate(mzO_p(0:np-1))
      allocate(mzNe_p(0:np-1))
      allocate(mzMg_p(0:np-1))
      allocate(mzSi_p(0:np-1))
      allocate(mzFe_p(0:np-1))
      allocate(mzZ_p(0:np-1))
      allocate(ts_p(0:np-1))
      allocate(flagc_p(0:np-1))
      allocate(flagfd_p(0:np-1))

end subroutine

subroutine allocate_dm(ndm)
      use gcdp_dm

      implicit none

      integer,intent(in) :: ndm

      if(allocated(id_dm)) then
        deallocate(id_dm)
        deallocate(nnb_dm)
        deallocate(m_dm)
        deallocate(x_dm)
        deallocate(y_dm)
        deallocate(z_dm)
        deallocate(vx_dm)
        deallocate(vy_dm)
        deallocate(vz_dm)
        deallocate(rho_dm)
        deallocate(h_dm)
        deallocate(list_adm)
      endif

      allocate(id_dm(0:ndm-1))
      allocate(nnb_dm(0:ndm-1))
      allocate(m_dm(0:ndm-1))
      allocate(x_dm(0:ndm-1))
      allocate(y_dm(0:ndm-1))
      allocate(z_dm(0:ndm-1))
      allocate(vx_dm(0:ndm-1))
      allocate(vy_dm(0:ndm-1))
      allocate(vz_dm(0:ndm-1))
      allocate(rho_dm(0:ndm-1))
      allocate(h_dm(0:ndm-1))
      allocate(list_adm(0:ndm-1))

end subroutine

! allocate general particle data
subroutine allocate_particle(np)
      use particle

      implicit none

      integer,intent(in) :: np

      if(allocated(idp)) then
        deallocate(idp)
        deallocate(listp)
        deallocate(xp)
        deallocate(yp)
        deallocate(massp)
        deallocate(metp)
        deallocate(hp)
      endif

      allocate(idp(0:np-1))
      allocate(listp(0:np-1))
      allocate(xp(0:np-1))
      allocate(yp(0:np-1))
      allocate(massp(0:np-1))
      allocate(metp(0:np-1))
      allocate(hp(0:np-1))
end subroutine

! allocate gtree
! reallocate gtree
subroutine allocate_gtree(ntr)
      use gcdp_gtree

      implicit none
      integer,intent(in) :: ntr

      if(allocated(np_gtr)) then
        deallocate(np_gtr)
        deallocate(pn_gtr)
        deallocate(l_gtr)
        deallocate(hm_gtr)
        deallocate(cx_gtr)
        deallocate(cy_gtr)
        deallocate(cz_gtr)
        deallocate(daughter_gtr)
        deallocate(next_gtr)
        deallocate(proc_gtr)
      endif

      allocate(np_gtr(0:ntr))
      allocate(pn_gtr(0:ntr))
      allocate(l_gtr(0:ntr))
      allocate(hm_gtr(0:ntr))
      allocate(cx_gtr(0:ntr))
      allocate(cy_gtr(0:ntr))
      allocate(cz_gtr(0:ntr))
      allocate(daughter_gtr(0:ntr))
      allocate(next_gtr(0:ntr))
      allocate(proc_gtr(0:ntr))

end subroutine
