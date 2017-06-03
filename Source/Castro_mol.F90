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

    use advection_util_module, only: compute_cfl, divu, normalize_species_fluxes, &
                                     ht, apply_av, construct_hydro_update, scale_flux
    use flatten_module, only: uflaten
    use riemann_module, only: cmpflx
    use ppm_module, only: ppm_reconstruct, ppm_int_profile
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

    integer :: ngf

    integer :: edge_lo(3), edge_hi(3)
    integer :: g_lo(3), g_hi(3)

    ngf = 1

    g_lo = lo - ngf
    g_hi = hi + ngf

    edge_lo = lo
    edge_hi = hi + 1

    ! Check if we have violated the CFL criterion.
    call compute_cfl(q, q_lo, q_hi, &
                     qaux, qa_lo, qa_hi, &
                     lo, hi, dt, dx, courno)

    ! Compute divergence of velocity field
    call divu(g_lo, g_hi, q, q_lo, q_hi, dx, h%div)

    ! Compute flattening coefficient for slope calculations.
    call uflaten(g_lo, g_hi, q, q_lo, q_hi, h)

    ! Create polynomial interpolation of fluid state.
    call ppm_reconstruct(q, q_lo, q_hi, h, g_lo, g_hi)

    ! Integrate under the reconstructed polynomial to get the edge state.
    call ppm_int_profile(h, g_lo, g_hi)

    ! Compute F^x
    call cmpflx(flux1, h%q1, f1_lo, f1_hi, &
                qaux, qa_lo, qa_hi, &
                h, 1, f1_lo, f1_hi, domlo, domhi)

    ! Compute F^y
    call cmpflx(flux2, h%q2, f2_lo, f2_hi, &
                qaux, qa_lo, qa_hi, &
                h, 2, f2_lo, f2_hi, domlo, domhi)

    ! Compute F^z
    call cmpflx(flux3, h%q3, f3_lo, f3_hi, &
                qaux, qa_lo, qa_hi, &
                h, 3, f3_lo, f3_hi, domlo, domhi)

    ! Apply artificial viscosity
    call apply_av(f1_lo, f1_hi, 1, dx, h, uin, uin_lo, uin_hi, flux1, f1_lo, f1_hi)
    call apply_av(f2_lo, f2_hi, 2, dx, h, uin, uin_lo, uin_hi, flux2, f2_lo, f2_hi)
    call apply_av(f3_lo, f3_hi, 3, dx, h, uin, uin_lo, uin_hi, flux3, f3_lo, f3_hi)

    ! Normalize species fluxes
    call normalize_species_fluxes(f1_lo, f1_hi, flux1, f1_lo, f1_hi)
    call normalize_species_fluxes(f2_lo, f2_hi, flux2, f2_lo, f2_hi)
    call normalize_species_fluxes(f3_lo, f3_hi, flux3, f3_lo, f3_hi)

    ! Scale the fluxes for the form we expect later in refluxing.
    call scale_flux(f1_lo, f1_hi, flux1, f1_lo, f1_hi, area1, a1_lo, a1_hi, dt)
    call scale_flux(f2_lo, f2_hi, flux2, f2_lo, f2_hi, area2, a2_lo, a2_hi, dt)
    call scale_flux(f3_lo, f3_hi, flux3, f3_lo, f3_hi, area3, a3_lo, a3_hi, dt)

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

end module mol_module
