module probdata_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), allocatable, public :: p_ambient, dens_ambient, exp_energy
  real(rt), allocatable, public :: r_init
  integer,  allocatable, public :: nsub
  integer,  allocatable, public :: probtype

#ifdef CUDA
  attributes(managed) :: p_ambient, dens_ambient, exp_energy, &
                         r_init, nsub, probtype
#endif
  
end module probdata_module
