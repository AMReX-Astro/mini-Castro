module network

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, parameter :: nspec = 13

  ! Number of nucleons per element
  real(rt), parameter :: aion(nspec) = [ 4.0d0, 12.0d0, 16.0d0, 20.0d0, 24.0d0, &
                                        28.0d0, 32.0d0, 36.0d0, 40.0d0, 44.0d0, &
                                        48.0d0, 52.0d0, 56.0d0]

  ! Number of protons per element
  real(rt), parameter :: zion(nspec) = [ 2.0d0,  6.0d0,  8.0d0, 10.0d0, 12.0d0, &
                                        14.0d0, 16.0d0, 18.0d0, 20.0d0, 22.0d0, &
                                        24.0d0, 26.0d0, 28.0d0]

  real(rt), parameter :: aion_inv(nspec) = 1.0d0 / aion(:)

end module network
