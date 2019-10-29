module reduction_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine reduce_add(x, y)

    !$acc routine seq

    implicit none

    ! Add y to x atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
    t = atomicAdd(x, y)
#else
    x = x + y
#endif

  end subroutine reduce_add



  CASTRO_FORT_DEVICE subroutine reduce_min(x, y)

    !$acc routine seq

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_ACC)
    t = atomicMin(x, y)
#else
    x = min(x, y)
#endif

  end subroutine reduce_min

end module reduction_module
