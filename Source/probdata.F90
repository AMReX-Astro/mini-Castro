module probdata_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), save ::  p_ambient, dens_ambient, exp_energy
  real(rt), save ::  r_init
  integer,  save ::  nsub
  integer,  save :: probtype, idir

end module probdata_module
