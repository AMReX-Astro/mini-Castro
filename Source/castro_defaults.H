
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

amrex::Real Castro::small_dens = -1.e200;
amrex::Real Castro::small_temp = -1.e200;
std::string Castro::xl_ext_bc_type = "";
std::string Castro::xr_ext_bc_type = "";
std::string Castro::yl_ext_bc_type = "";
std::string Castro::yr_ext_bc_type = "";
std::string Castro::zl_ext_bc_type = "";
std::string Castro::zr_ext_bc_type = "";
amrex::Real Castro::fixed_dt = -1.0;
amrex::Real Castro::cfl = 0.8;
int         Castro::do_acc = -1;
#ifdef DEBUG
int         Castro::print_fortran_warnings = 1;
#else
int         Castro::print_fortran_warnings = 0;
#endif
int         Castro::sum_interval = -1;
amrex::Real Castro::sum_per = -1.0e0;
int         Castro::hard_cfl_limit = 1;
