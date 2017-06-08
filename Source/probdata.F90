module probdata_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), allocatable :: p_ambient, dens_ambient, exp_energy
  real(rt), allocatable :: r_init
  integer,  allocatable :: nsub
  integer,  allocatable :: probtype, idir

#ifdef CUDA
  attributes(managed) :: p_ambient, dens_ambient, exp_energy, r_init, nsub, probtype, idir
#endif

end module probdata_module
