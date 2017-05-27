
module actual_network

  use bl_types

  implicit none

  integer, parameter :: nspec = 1
  integer, parameter :: nspec_evolve = 1
  integer, parameter :: naux =  0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len=16), save :: aux_names(naux)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save   :: aion(nspec), zion(nspec)

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 1

contains
  
  subroutine actual_network_init

    spec_names(1) = "X"

    short_spec_names(1) = "X"

    aion(1) = 1.0

    zion(1) = 1.0

  end subroutine actual_network_init

end module actual_network
