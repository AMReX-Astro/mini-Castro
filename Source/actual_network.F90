
module actual_network

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: nspec = 1
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), zion(:)

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 1

#ifdef CUDA
  attributes(managed) :: aion, zion
#endif

contains
  
  subroutine actual_network_init

    spec_names(1) = "X"

    short_spec_names(1) = "X"

    allocate(aion(nspec))
    allocate(zion(nspec))

    aion(1) = 1.0

    zion(1) = 1.0

  end subroutine actual_network_init

end module actual_network
