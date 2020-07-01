module reduction_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  CASTRO_FORT_DEVICE subroutine reduce_add(x, y)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    ! Add y to x atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_ACC) && !defined(AMREX_USE_OMP_OFFLOAD)
    t = atomicAdd(x, y)
#else
#ifdef AMREX_USE_OMP
    !$omp atomic
#endif
    x = x + y
#endif

  end subroutine reduce_add



  CASTRO_FORT_DEVICE subroutine reduce_min(x, y)

#ifdef AMREX_USE_ACC
    !$acc routine seq
#endif

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(rt), intent(in   ) :: y
    real(rt), intent(inout) :: x

    real(rt) :: t

#ifdef AMREX_USE_OMP_OFFLOAD
    !$omp declare target
#endif

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_ACC) && !defined(AMREX_USE_OMP_OFFLOAD)
    t = atomicMin(x, y)
#else
#ifdef AMREX_USE_OMP
    !$omp atomic
#endif
    x = min(x, y)
#endif

  end subroutine reduce_min

end module reduction_module
