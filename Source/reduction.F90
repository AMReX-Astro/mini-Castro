module reduction_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  AMREX_CUDA_FORT_DEVICE subroutine reduce_add(x, y)

    implicit none

    ! Add y to x atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#ifdef AMREX_USE_CUDA
    t = atomicAdd(x, y)
#else
    x = x + y
#endif

  end subroutine reduce_add



  AMREX_CUDA_FORT_DEVICE subroutine reduce_min(x, y)

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#ifdef AMREX_USE_CUDA
    t = atomicMin(x, y)
#else
    x = min(x, y)
#endif

  end subroutine reduce_min

end module reduction_module
