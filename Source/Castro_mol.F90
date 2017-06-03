! advection routines in support of method of lines integration

module mol_module

  implicit none

contains
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine mol_single_stage(time, &
                              lo, hi, domlo, domhi, &
                              uin, uin_lo, uin_hi, &
                              uout, uout_lo, uout_hi, &
                              q, q_lo, q_hi, &
                              qaux, qa_lo, qa_hi, &
                              update, updt_lo, updt_hi, &
                              dx, dt, h, &
                              flux1, f1_lo, f1_hi, &
                              flux2, f2_lo, f2_hi, &
                              flux3, f3_lo, f3_hi, &
                              area1, a1_lo, a1_hi, &
                              area2, a2_lo, a2_hi, &
                              area3, a3_lo, a3_hi, &
                              vol, vol_lo, vol_hi, &
                              courno, verbose)

    use advection_util_module, only: ht, construct_hydro_update
    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, QVAR, NVAR, NQAUX

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), verbose
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: updt_lo(3), updt_hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)
    integer,  intent(in   ) :: a1_lo(3), a1_hi(3)
    integer,  intent(in   ) :: a2_lo(3), a2_hi(3)
    integer,  intent(in   ) :: a3_lo(3), a3_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(inout) :: flux1(f1_lo(1):f1_hi(1), f1_lo(2):f1_hi(2), f1_lo(3):f1_hi(3), NVAR)
    real(rt), intent(inout) :: flux2(f2_lo(1):f2_hi(1), f2_lo(2):f2_hi(2), f2_lo(3):f2_hi(3), NVAR)
    real(rt), intent(inout) :: flux3(f3_lo(1):f3_hi(1), f3_lo(2):f3_hi(2), f3_lo(3):f3_hi(3), NVAR)
    real(rt), intent(in   ) :: area1(a1_lo(1):a1_hi(1), a1_lo(2):a1_hi(2), a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: area2(a2_lo(1):a2_hi(1), a2_lo(2):a2_hi(2), a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: area3(a3_lo(1):a3_hi(1), a3_lo(2):a3_hi(2), a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(in   ) :: dx(3), dt, time
    real(rt), intent(inout) :: courno

    type(ht), intent(inout) :: h

    integer, parameter :: ngf = 1

    integer :: g_lo(3), g_hi(3)

    g_lo = lo - ngf
    g_hi = hi + ngf

    ! Construct edge states as inputs to flux construction

    call prepare_for_fluxes(g_lo, g_hi, dt, dx, courno, h, &
                            q, q_lo, q_hi, &
                            qaux, qa_lo, qa_hi)

    ! Compute F^x

    call construct_flux(f1_lo, f1_hi, domlo, domhi, h, dx, dt, 1, &
                        uin, uin_lo, uin_hi, &
                        flux1, h%q1, f1_lo, f1_hi, &
                        area1, a1_lo, a1_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^y

    call construct_flux(f2_lo, f2_hi, domlo, domhi, h, dx, dt, 2, &
                        uin, uin_lo, uin_hi, &
                        flux2, h%q2, f2_lo, f2_hi, &
                        area2, a2_lo, a2_hi, &
                        qaux, qa_lo, qa_hi)

    ! Compute F^z

    call construct_flux(f3_lo, f3_hi, domlo, domhi, h, dx, dt, 3, &
                        uin, uin_lo, uin_hi, &
                        flux3, h%q3, f3_lo, f3_hi, &
                        area3, a3_lo, a3_hi, &
                        qaux, qa_lo, qa_hi)

    ! Create an update source term based on the flux divergence.
    call construct_hydro_update(lo, hi, dx, dt, h, &
                                flux1, f1_lo, f1_hi, &
                                flux2, f2_lo, f2_hi, &
                                flux3, f3_lo, f3_hi, &
                                area1, a1_lo, a1_hi, &
                                area2, a2_lo, a2_hi, &
                                area3, a3_lo, a3_hi, &
                                vol, vol_lo, vol_hi, &
                                update, updt_lo, updt_hi)

  end subroutine mol_single_stage


  subroutine prepare_for_fluxes(lo, hi, dt, dx, courno, h, &
                                q, q_lo, q_hi, &
                                qaux, qa_lo, qa_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NQ, NQAUX
    use advection_util_module, only: ht, compute_cfl, divu
    use flatten_module, only: uflaten
    use ppm_module, only: ppm_reconstruct, ppm_int_profile

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)

    real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: dx(3), dt
    real(rt), intent(inout) :: courno

    type(ht), intent(inout) :: h

    ! Check if we have violated the CFL criterion.
    call compute_cfl(lo, hi, dt, dx, courno, &
                     q, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi)

    ! Compute divergence of velocity field.
    call divu(lo, hi, dx, q, q_lo, q_hi, h%div)

    ! Compute flattening coefficient for slope calculations.
    call uflaten(lo, hi, q, q_lo, q_hi, h)

    ! Create polynomial interpolation of fluid state.
    call ppm_reconstruct(lo, hi, q, q_lo, q_hi, h)

    ! Integrate under the reconstructed polynomial to get the edge state.
    call ppm_int_profile(lo, hi, h)

  end subroutine prepare_for_fluxes



  subroutine construct_flux(lo, hi, domlo, domhi, h, dx, dt, idir, &
                            uin, uin_lo, uin_hi, &
                            flux, qint, f_lo, f_hi, &
                            area, a_lo, a_hi, &
                            qaux, qa_lo, qa_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX
    use advection_util_module, only: ht, apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3), idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(inout) :: qint(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: dx(3), dt

    type(ht), intent(inout) :: h

    call cmpflx(lo, hi, domlo, domhi, h, idir, flux, qint, f_lo, f_hi, qaux, qa_lo, qa_hi)
    call apply_av(lo, hi, idir, dx, h, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine construct_flux

end module mol_module
