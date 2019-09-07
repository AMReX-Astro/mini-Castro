module eos_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public eos_init, eos

contains

  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init() bind(c, name='eos_init')

    use actual_eos_module, only: actual_eos_init

    implicit none

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

    !$gpu

    ! Get abar, zbar, etc.

    call composition(state)

    ! Call the EOS.

    call actual_eos(input, state)

  end subroutine eos



  subroutine eos_finalize() bind(c, name='eos_finalize')

    use actual_eos_module, only: actual_eos_finalize

    implicit none

    call actual_eos_finalize()

  end subroutine eos_finalize

end module eos_module
