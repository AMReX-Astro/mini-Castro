
! This module stores the runtime parameters that define the problem domain.  
! These parameter are initialized in set_problem_params().

module prob_params_module

  use meth_params_module, only: UMX, UMZ
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  ! geometry information
  real(rt), allocatable :: problo(:), probhi(:)

  ! indices that we use for dimension agnostic routines 
  ! to ensure we don't illegally access non-existent ghost cells
  integer, allocatable :: dg(:)

  ! grid information
  integer         , save        :: max_level
  integer         , allocatable :: domlo_level(:,:)
  integer         , allocatable :: domhi_level(:,:)
  integer         , allocatable :: ref_ratio(:,:)
  integer         , allocatable :: n_error_buf(:)
  integer         , allocatable :: blocking_factor(:)

  integer, parameter :: MAX_MOM_INDEX = 5

  type momflux_t
     ! we want this to be able to use UMX, UMY, and UMZ to index here, but
     ! we can't use those to allocate, since they are not know until runtime.
     ! dynamic allocation might mess with GPUs, so we make this big enough
     ! to definitely contain UMX, UMY, and UMZ, and then check this when
     ! we fill it
     logical :: comp(MAX_MOM_INDEX)
  end type momflux_t

  ! one component for each coordinate direction flux
  type (momflux_t), save :: mom_flux_has_p(3)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: dg
  attributes(managed) :: problo, probhi
  attributes(managed) :: domlo_level, domhi_level
#endif
  
end module prob_params_module
