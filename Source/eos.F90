module eos_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public eos_init, eos

contains

  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init() bind(c, name='eos_init')

    use eos_type_module, only: mintemp, maxtemp, mindens, maxdens, minx, maxx, &
                               minye, maxye, mine, maxe, minp, maxp
    use actual_eos_module, only: actual_eos_init

    implicit none

    ! Allocate and set default values

    allocate(mintemp)
    allocate(maxtemp)
    allocate(mindens)
    allocate(maxdens)
    allocate(minx)
    allocate(maxx)
    allocate(minye)
    allocate(maxye)
    allocate(mine)
    allocate(maxe)
    allocate(minp)
    allocate(maxp)

    mintemp = 1.d-200
    maxtemp = 1.d200
    mindens = 1.d-200
    maxdens = 1.d200
    minx    = 1.d-200
    maxx    = 1.d0 + 1.d-12
    minye   = 1.d-200
    maxye   = 1.d0 + 1.d-12
    mine    = 1.d-200
    maxe    = 1.d200
    minp    = 1.d-200
    maxp    = 1.d200

    ! Set up any specific parameters or initialization steps required by the EOS we are using.

    call actual_eos_init()

  end subroutine eos_init



  subroutine eos(input, state)

    use eos_type_module, only: eos_t, composition
    use actual_eos_module, only: actual_eos

    implicit none

    ! Input arguments

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    logical :: has_been_reset

    !$gpu

    ! Get abar, zbar, etc.

    call composition(state)

    ! Force the inputs to be valid.

    has_been_reset = .false.
    call reset_inputs(input, state, has_been_reset)

    ! Call the EOS.

    if (.not. has_been_reset) then
       call actual_eos(input, state)
    endif

  end subroutine eos



  subroutine reset_inputs(input, state, has_been_reset)

    use eos_type_module, only: eos_t, eos_input_rt, eos_input_re, eos_input_rp

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    ! Reset the input quantities to valid values. For inputs other than rho and T,
    ! this will evolve an EOS call, which will negate the need to do the main EOS call.

    if (input .eq. eos_input_rt) then

       call reset_rho(state, has_been_reset)
       call reset_T(state, has_been_reset)

    elseif (input .eq. eos_input_rp) then

       call reset_rho(state, has_been_reset)
       call reset_p(state, has_been_reset)

    elseif (input .eq. eos_input_re) then

       call reset_rho(state, has_been_reset)
       call reset_e(state, has_been_reset)

    endif

  end subroutine reset_inputs



  ! For density, just ensure that it is within mindens and maxdens.

  subroutine reset_rho(state, has_been_reset)

    use eos_type_module, only: eos_t, mindens, maxdens

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine reset_rho



  ! For temperature, just ensure that it is within mintemp and maxtemp.

  subroutine reset_T(state, has_been_reset)

    use eos_type_module, only: eos_t, mintemp, maxtemp

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % T = min(maxtemp, max(mintemp, state % T))

  end subroutine reset_T



  subroutine reset_e(state, has_been_reset)

    use eos_type_module, only: eos_t, mine, maxe

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % e .lt. mine .or. state % e .gt. maxe) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_e


  subroutine reset_p(state, has_been_reset)

    use eos_type_module, only: eos_t, minp, maxp

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % p .lt. minp .or. state % p .gt. maxp) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_p



  ! Given an EOS state, ensure that rho and T are
  ! valid, then call with eos_input_rt.

  subroutine eos_reset(state, has_been_reset)

    use actual_eos_module, only: actual_eos
    use eos_type_module, only: eos_t, eos_input_rt, mintemp, maxtemp, mindens, maxdens

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

    call actual_eos(eos_input_rt, state)

    has_been_reset = .true.

  end subroutine eos_reset



  subroutine eos_finalize() bind(c, name='eos_finalize')

    use eos_type_module, only: mintemp, maxtemp, mindens, maxdens, &
                               minx, maxx, minye, maxye, &
                               mine, maxe, minp, maxp, &
                               mins, maxs, minh, maxh
    use actual_eos_module, only: actual_eos_finalize

    implicit none

    deallocate(mintemp)
    deallocate(maxtemp)
    deallocate(mindens)
    deallocate(maxdens)
    deallocate(minx)
    deallocate(maxx)
    deallocate(minye)
    deallocate(maxye)
    deallocate(mine)
    deallocate(maxe)
    deallocate(minp)
    deallocate(maxp)

    call actual_eos_finalize()

  end subroutine eos_finalize

end module eos_module
