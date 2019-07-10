! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
!
! aion_inv   -- 1/aion
!
! spec_names -- the name of the isotope
! short_spec_names -- the abbreviated name of the isotope
!
! This module contains the following routines:
!
!  network_init          -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name

module network

  use amrex_fort_module, only: rt => amrex_real
  use actual_network

  implicit none

  logical :: network_initialized = .false.

  ! this will be computed here, not in the actual network
  real(rt), allocatable :: aion_inv(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion_inv
#endif

contains

  subroutine network_init() bind(c, name='network_init')

    use amrex_error_module, only: amrex_error
    use amrex_constants_module, only: ONE

    implicit none

    allocate(aion_inv(nspec))

    ! Call the specific network initialization.

    call actual_network_init()

    aion_inv(:) = ONE / aion(:)

    network_initialized = .true.

  end subroutine network_init


  subroutine network_finalize() bind(c, name='network_finalize')

    use actual_network, only: actual_network_finalize

    implicit none

    deallocate(aion_inv)

    call actual_network_finalize()

  end subroutine network_finalize

end module network
