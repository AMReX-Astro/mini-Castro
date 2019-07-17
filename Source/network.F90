module network

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: nspec = 13

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: is32  = 7
  integer, parameter :: iar36 = 8
  integer, parameter :: ica40 = 9
  integer, parameter :: iti44 = 10
  integer, parameter :: icr48 = 11
  integer, parameter :: ife52 = 12
  integer, parameter :: ini56 = 13

  real(rt), allocatable :: aion(:), zion(:), aion_inv(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, aion_inv
#endif

  character(len=5), save :: short_spec_names(nspec)

contains

  subroutine network_init() bind(c, name='network_init')

    implicit none

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(aion_inv(nspec))

    ! Set the number of nucleons in the element
    aion(ihe4)  = 4.0d0
    aion(ic12)  = 12.0d0
    aion(io16)  = 16.0d0
    aion(ine20) = 20.0d0
    aion(img24) = 24.0d0
    aion(isi28) = 28.0d0
    aion(is32)  = 32.0d0
    aion(iar36) = 36.0d0
    aion(ica40) = 40.0d0
    aion(iti44) = 44.0d0
    aion(icr48) = 48.0d0
    aion(ife52) = 52.0d0
    aion(ini56) = 56.0d0

    ! Set the number of protons in the element
    zion(ihe4)  = 2.0d0
    zion(ic12)  = 6.0d0
    zion(io16)  = 8.0d0
    zion(ine20) = 10.0d0
    zion(img24) = 12.0d0
    zion(isi28) = 14.0d0
    zion(is32)  = 16.0d0
    zion(iar36) = 18.0d0
    zion(ica40) = 20.0d0
    zion(iti44) = 22.0d0
    zion(icr48) = 24.0d0
    zion(ife52) = 26.0d0
    zion(ini56) = 28.0d0

    aion_inv(:) = 1.d0 / aion(:)

    ! Set the species names
    short_spec_names(ihe4)  = 'he4'
    short_spec_names(ic12)  = 'c12'
    short_spec_names(io16)  = 'o16'
    short_spec_names(ine20) = 'ne20'
    short_spec_names(img24) = 'mg24'
    short_spec_names(isi28) = 'si28'
    short_spec_names(is32)  = 's32'
    short_spec_names(iar36) = 'ar36'
    short_spec_names(ica40) = 'ca40'
    short_spec_names(iti44) = 'ti44'
    short_spec_names(icr48) = 'cr48'
    short_spec_names(ife52) = 'fe52'
    short_spec_names(ini56) = 'ni56'

  end subroutine network_init


  subroutine network_finalize() bind(c, name='network_finalize')

    implicit none

    deallocate(aion)
    deallocate(zion)
    deallocate(aion_inv)

  end subroutine network_finalize

end module network
