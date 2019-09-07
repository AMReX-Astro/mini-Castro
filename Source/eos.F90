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

    !$gpu

    ! Get abar, zbar, etc.

    call composition(state)

    ! Call the EOS.

    call actual_eos(input, state)

  end subroutine eos



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
