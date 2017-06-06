! advection routines in support of method of lines integration

module mol_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  subroutine prepare_for_fluxes(lo, hi, dt, dx, courno, &
                                q, flatn, q_lo, q_hi, &
                                div, g_lo, g_hi, &
                                qaux, qa_lo, qa_hi, &
                                sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                                qm, qp, It_lo, It_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use advection_util_module, only: compute_cfl, divu
    use flatten_module, only: uflaten
    use ppm_module, only: ppm_reconstruct

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: g_lo(3), g_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: st_lo(3), st_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: flatn(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3))
    real(rt), intent(inout) :: div(g_lo(1):g_hi(1),g_lo(2):g_hi(2),g_lo(3):g_hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: sxm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sxp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sym(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: syp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(in   ) :: dx(3), dt
    real(rt), intent(inout) :: courno

    ! Check if we have violated the CFL criterion.
    call compute_cfl(lo, hi, dt, dx, courno, &
                     q, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi)

    ! Compute divergence of velocity field.
    call divu(lo, hi, dx, q, q_lo, q_hi, div, g_lo, g_hi)

    ! Compute flattening coefficient for slope calculations.
    call uflaten(lo, hi, q, flatn, q_lo, q_hi)

    ! Create polynomial interpolation of fluid state.
    call ppm_reconstruct(lo, hi, q, flatn, q_lo, q_hi, &
                         sxm, sxp, sym, syp, szm, szp, st_lo, st_hi)

  end subroutine prepare_for_fluxes



#ifdef CUDA
  attributes(device) &
#endif
  subroutine prepare_profile(lo, hi, &
                             q, flatn, q_lo, q_hi, &
                             sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                             qm, qp, It_lo, It_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use ppm_module, only: ppm_int_profile

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: st_lo(3), st_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: flatn(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3))
    real(rt), intent(inout) :: sxm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sxp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: sym(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: syp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szm(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: szp(st_lo(1):st_hi(1),st_lo(2):st_hi(2),st_lo(3):st_hi(3),NQ)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)

    ! Integrate under the reconstructed polynomial to get the edge state.
    call ppm_int_profile(lo, hi, sxm, sxp, sym, syp, szm, szp, st_lo, st_hi, &
                         qm, qp, It_lo, It_hi)

  end subroutine prepare_profile



#ifdef CUDA
  attributes(device) &
#endif
  subroutine construct_flux(lo, hi, domlo, domhi, dx, dt, idir, &
                            div, g_lo, g_hi, &
                            uin, uin_lo, uin_hi, &
                            qm, qp, It_lo, It_hi, &
                            flux, qint, f_lo, f_hi, &
                            area, a_lo, a_hi, &
                            qaux, qa_lo, qa_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use advection_util_module, only: apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: g_lo(3), g_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: It_lo(3), It_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: div(g_lo(1):g_hi(1), g_lo(2):g_hi(2), g_lo(3):g_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: qm(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(It_lo(1):It_hi(1),It_lo(2):It_hi(2),It_lo(3):It_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: dx(3), dt

    call cmpflx(lo, hi, domlo, domhi, idir, qm, qp, It_lo, It_hi, flux, qint, f_lo, f_hi, qaux, qa_lo, qa_hi)
    call apply_av(lo, hi, idir, dx, div, g_lo, g_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine construct_flux

end module mol_module
