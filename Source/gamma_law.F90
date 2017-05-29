! This is a constant gamma equation of state, using an ideal gas.

module actual_eos_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  character (len=64) :: eos_name = "gamma_law"

  double precision, save :: gamma_const

  logical, save :: assume_neutral

  ! boltzmann's constant
  real(rt), parameter :: k_B = 1.3806488e-16_rt   ! erg/K

  ! avogradro's Number
  real(rt), parameter :: n_A = 6.02214129e23_rt   ! mol^-1

#ifdef CUDA
  double precision, device :: gamma_const_d
  logical, device :: assume_neutral_d
#endif

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_gamma, eos_assume_neutral
    use bl_constants_module, only: ZERO

    implicit none

    ! constant ratio of specific heats
    if (eos_gamma .gt. ZERO) then
       gamma_const = eos_gamma
    else
       call bl_error("gamma_const cannot be < 0")
    end if

    assume_neutral = eos_assume_neutral

#ifdef CUDA
    gamma_const_d = gamma_const
    assume_neutral_d = assume_neutral
#endif

  end subroutine actual_eos_init



#ifdef CUDA
  attributes(device) &
#endif
  subroutine actual_eos(input, state)

    use eos_type_module, only: eos_t, &
                               eos_input_rt, eos_input_re, eos_input_rh, eos_input_rp, &
                               eos_input_th, eos_input_tp, eos_input_ps, eos_input_ph
    use bl_constants_module, only: ZERO, ONE
#if !(defined(ACC) || defined(CUDA))
    use bl_error_module, only: bl_error
#endif
#ifdef CUDA
    use network, only: aion => aion_d, zion => zion_d
#else
    use network, only: aion, zion
#endif

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    double precision, parameter :: R = k_B*n_A

    double precision :: poverrho

    double precision :: gamma
    logical :: neutral

#ifdef CUDA
    gamma = gamma_const_d
    neutral = assume_neutral_d
#else
    gamma = gamma_const
    neutral = assume_neutral
#endif

    ! Calculate mu.

    if (neutral) then
       state % mu = state % abar
    else
       state % mu = ONE / sum( (ONE + zion(:)) * state % xn(:) / aion(:) )
    endif

    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs
       state % cv = R / (state % mu * (gamma-ONE))
       state % e = state % cv * state % T
       state % p = (gamma-ONE) * state % rho * state % e
       state % gam1 = gamma

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: eos_input_rh is not supported in this EOS.')
#endif

    case (eos_input_tp)

       ! temp, pres, and xmass are inputs

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: eos_input_tp is not supported in this EOS.')
#endif

    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       poverrho = state % p / state % rho
       state % T = poverrho * state % mu * (ONE/R)
       state % e = poverrho * (ONE/(gamma-ONE))
       state % gam1 = gamma

    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       poverrho = (gamma - ONE) * state % e

       state % p = poverrho * state % rho
       state % T = poverrho * state % mu * (ONE/R)
       state % gam1 = gamma

       ! sound speed
       state % cs = sqrt(gamma * poverrho)

       state % dpdr_e = poverrho
       state % dpde = (gamma-ONE) * state % rho

       ! Try to avoid the expensive log function.  Since we don't need entropy
       ! in hydro solver, set it to an invalid but "nice" value for the plotfile.
       state % s = ONE

    case (eos_input_ps)

       ! pressure entropy, and xmass are inputs

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: eos_input_ps is not supported in this EOS.')
#endif

    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: eos_input_ph is not supported in this EOS.')
#endif

    case (eos_input_th)

       ! temperature, enthalpy and xmass are inputs

       ! This system is underconstrained.

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: eos_input_th is not a valid input for the gamma law EOS.')
#endif

    case default

#if !(defined(ACC) || defined(CUDA))
       call bl_error('EOS: invalid input.')
#endif

    end select

    ! Give dpdr a value for the purposes of the composition_derivatives routine.

    state % dPdr = ZERO

  end subroutine actual_eos

end module actual_eos_module
