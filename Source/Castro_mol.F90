! advection routines in support of method of lines integration

module mol_module

  use amrex_fort_module, only: rt => amrex_real, get_loop_bounds

  implicit none

contains

  AMREX_LAUNCH subroutine ca_prepare_for_fluxes(lo, hi, dx, dt, &
                                                q, q_lo, q_hi, &
                                                qaux, qa_lo, qa_hi, &
                                                flatn, f_lo, f_hi, &
                                                div, d_lo, d_hi, &
                                                sxm, sxm_lo, sxm_hi, &
                                                sxp, sxp_lo, sxp_hi, &
                                                sym, sym_lo, sym_hi, &
                                                syp, syp_lo, syp_hi, &
                                                szm, szm_lo, szm_hi, &
                                                szp, szp_lo, szp_hi) &
                                                bind(c,name='ca_prepare_for_fluxes')

    use meth_params_module, only: NQ, NQAUX
    use advection_util_module, only: compute_cfl, divu
    use flatten_module, only: uflaten
    use ppm_module, only: ppm_reconstruct

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    integer,  intent(in   ) :: sxm_lo(3), sxm_hi(3)
    integer,  intent(in   ) :: sxp_lo(3), sxp_hi(3)
    integer,  intent(in   ) :: sym_lo(3), sym_hi(3)
    integer,  intent(in   ) :: syp_lo(3), syp_hi(3)
    integer,  intent(in   ) :: szm_lo(3), szm_hi(3)
    integer,  intent(in   ) :: szp_lo(3), szp_hi(3)

    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    real(rt), intent(inout) :: div(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))
    real(rt), intent(inout) :: sxm(sxm_lo(1):sxm_hi(1),sxm_lo(2):sxm_hi(2),sxm_lo(3):sxm_hi(3),NQ)
    real(rt), intent(inout) :: sxp(sxp_lo(1):sxp_hi(1),sxp_lo(2):sxp_hi(2),sxp_lo(3):sxp_hi(3),NQ)
    real(rt), intent(inout) :: sym(sym_lo(1):sym_hi(1),sym_lo(2):sym_hi(2),sym_lo(3):sym_hi(3),NQ)
    real(rt), intent(inout) :: syp(syp_lo(1):syp_hi(1),syp_lo(2):syp_hi(2),syp_lo(3):syp_hi(3),NQ)
    real(rt), intent(inout) :: szm(szm_lo(1):szm_hi(1),szm_lo(2):szm_hi(2),szm_lo(3):szm_hi(3),NQ)
    real(rt), intent(inout) :: szp(szp_lo(1):szp_hi(1),szp_lo(2):szp_hi(2),szp_lo(3):szp_hi(3),NQ)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    ! Compute divergence of velocity field.
    call divu(blo, bhi, dx, q, q_lo, q_hi, div, d_lo, d_hi)

    ! Compute flattening coefficient for slope calculations.
    call uflaten(blo, bhi, q, q_lo, q_hi, flatn, f_lo, f_hi)

    ! Create polynomial interpolation of fluid state.
    call ppm_reconstruct(blo, bhi, &
                         q, q_lo, q_hi, &
                         flatn, f_lo, f_hi, &
                         sxm, sxm_lo, sxm_hi, &
                         sxp, sxp_lo, sxp_hi, &
                         sym, sym_lo, sym_hi, &
                         syp, syp_lo, syp_hi, &
                         szm, szm_lo, szm_hi, &
                         szp, szp_lo, szp_hi)

  end subroutine ca_prepare_for_fluxes



  AMREX_LAUNCH subroutine ca_prepare_profile(lo, hi, &
                                             q, q_lo, q_hi, &
                                             sxm, sxm_lo, sxm_hi, &
                                             sxp, sxp_lo, sxp_hi, &
                                             sym, sym_lo, sym_hi, &
                                             syp, syp_lo, syp_hi, &
                                             szm, szm_lo, szm_hi, &
                                             szp, szp_lo, szp_hi, &
                                             qm, qm_lo, qm_hi, &
                                             qp, qp_lo, qp_hi) &
                                             bind(c,name='ca_prepare_profile')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use ppm_module, only: ppm_int_profile

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: sxm_lo(3), sxm_hi(3)
    integer,  intent(in   ) :: sxp_lo(3), sxp_hi(3)
    integer,  intent(in   ) :: sym_lo(3), sym_hi(3)
    integer,  intent(in   ) :: syp_lo(3), syp_hi(3)
    integer,  intent(in   ) :: szm_lo(3), szm_hi(3)
    integer,  intent(in   ) :: szp_lo(3), szp_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)

    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(in   ) :: sxm(sxm_lo(1):sxm_hi(1),sxm_lo(2):sxm_hi(2),sxm_lo(3):sxm_hi(3),NQ)
    real(rt), intent(in   ) :: sxp(sxp_lo(1):sxp_hi(1),sxp_lo(2):sxp_hi(2),sxp_lo(3):sxp_hi(3),NQ)
    real(rt), intent(in   ) :: sym(sym_lo(1):sym_hi(1),sym_lo(2):sym_hi(2),sym_lo(3):sym_hi(3),NQ)
    real(rt), intent(in   ) :: syp(syp_lo(1):syp_hi(1),syp_lo(2):syp_hi(2),syp_lo(3):syp_hi(3),NQ)
    real(rt), intent(in   ) :: szm(szm_lo(1):szm_hi(1),szm_lo(2):szm_hi(2),szm_lo(3):szm_hi(3),NQ)
    real(rt), intent(in   ) :: szp(szp_lo(1):szp_hi(1),szp_lo(2):szp_hi(2),szp_lo(3):szp_hi(3),NQ)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    ! Integrate under the reconstructed polynomial to get the edge state.
    call ppm_int_profile(blo, bhi, &
                         sxm, sxm_lo, sxm_hi, &
                         sxp, sxp_lo, sxp_hi, &
                         sym, sym_lo, sym_hi, &
                         syp, syp_lo, syp_hi, &
                         szm, szm_lo, szm_hi, &
                         szp, szp_lo, szp_hi, &
                         qm, qm_lo, qm_hi, &
                         qp, qp_lo, qp_hi)

  end subroutine ca_prepare_profile



  AMREX_LAUNCH subroutine ca_construct_flux(lo, hi, domlo, domhi, dx, dt, idir, &
                                            uin, uin_lo, uin_hi, &
                                            div, div_lo, div_hi, &
                                            qaux, qa_lo, qa_hi, &
                                            qm, qm_lo, qm_hi, &
                                            qp, qp_lo, qp_hi, &
                                            qint, qe_lo, qe_hi, &
                                            flux, f_lo, f_hi, &
                                            area, a_lo, a_hi) &
                                            bind(c,name='ca_construct_flux')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use advection_util_module, only: apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(in   ) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    call cmpflx(blo, bhi, domlo, domhi, idir, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                qint, qe_lo, qe_hi, flux, f_lo, f_hi, qaux, qa_lo, qa_hi)
    call apply_av(blo, bhi, idir, dx, div, div_lo, div_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(blo, bhi, flux, f_lo, f_hi)
    call scale_flux(blo, bhi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine ca_construct_flux

end module mol_module
