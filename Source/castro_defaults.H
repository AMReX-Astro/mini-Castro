
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

int         Castro::state_interp_order = 1;
int         Castro::lin_limit_state_interp = 0;
int         Castro::state_nghost = 0;
amrex::Real Castro::small_dens = -1.e200;
amrex::Real Castro::small_temp = -1.e200;
int         Castro::ppm_predict_gammae = 0;
int         Castro::ppm_reference_eigenvectors = 0;
int         Castro::hybrid_riemann = 0;
int         Castro::use_colglaz = -1;
int         Castro::riemann_solver = 0;
int         Castro::cg_maxiter = 12;
amrex::Real Castro::cg_tol = 1.0e-5;
int         Castro::cg_blend = 2;
int         Castro::use_flattening = 1;
int         Castro::dual_energy_update_E_from_e = 1;
amrex::Real Castro::dual_energy_eta1 = 1.0e0;
amrex::Real Castro::dual_energy_eta2 = 1.0e-4;
amrex::Real Castro::dual_energy_eta3 = 1.0e0;
int         Castro::limit_fluxes_on_small_dens = 0;
int         Castro::density_reset_method = 1;
int         Castro::allow_negative_energy = 0;
int         Castro::allow_small_energy = 1;
int         Castro::first_order_hydro = 0;
std::string Castro::xl_ext_bc_type = "";
std::string Castro::xr_ext_bc_type = "";
std::string Castro::yl_ext_bc_type = "";
std::string Castro::yr_ext_bc_type = "";
std::string Castro::zl_ext_bc_type = "";
std::string Castro::zr_ext_bc_type = "";
int         Castro::hse_zero_vels = 0;
int         Castro::hse_interp_temp = 0;
int         Castro::hse_reflect_vels = 0;
amrex::Real Castro::fixed_dt = -1.0;
amrex::Real Castro::initial_dt = -1.0;
amrex::Real Castro::dt_cutoff = 0.0;
amrex::Real Castro::max_dt = 1.e200;
amrex::Real Castro::cfl = 0.8;
amrex::Real Castro::init_shrink = 1.0;
amrex::Real Castro::change_max = 1.1;
int         Castro::do_acc = -1;
#ifdef DEBUG
int         Castro::print_fortran_warnings = 1;
#else
int         Castro::print_fortran_warnings = 0;
#endif
int         Castro::sum_interval = -1;
amrex::Real Castro::sum_per = -1.0e0;
int         Castro::hard_cfl_limit = 1;
std::string Castro::job_name = "";
