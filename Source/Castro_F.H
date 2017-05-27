#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_network_init();

  void burner_init();

  void ca_extern_init(const int* name, const int* namlen);

  void ca_get_num_spec(int* nspec);

  void ca_get_num_aux(int* naux);

  void ca_get_spec_names(int* spec_names, int* ispec, int* len);

  void ca_get_spec_az(int* ispec, amrex::Real* A, amrex::Real* Z);

  void ca_get_aux_names(int* aux_names, int* iaux, int* len);

  void ca_get_qvar(int* qvar);

  void ca_get_nqaux(int* nqaux);

  void ca_set_special_tagging_flag(const amrex::Real& max_den, const int* flag);

  void ca_set_amr_info(const int& level, const int& iteration, const int& ncycle,
		       const amrex::Real& time, const amrex::Real& dt);

  void ca_get_method_params(int* HYP_GROW);

  void ca_set_method_params
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv);


  void ca_set_castro_method_params();

  void ca_init_godunov_indices();
  void ca_set_problem_params
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Interior_value, const int& Inflow_value,
     const int& Outflow_value, const int& Symmetry_value,
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const amrex::Real* problo, const amrex::Real* probhi,
     const amrex::Real* center);

  void ca_set_grid_info
    (const int& max_level, const amrex::Real* dx_level,
     const int* domlo_level, const int* domhi_level,
     const int* ref_ratio, const int* n_error_buf,
     const int* blocking_factor);

  void ca_get_tagging_params(const int* name, const int* namlen);

  void ca_ctoprim(const int* lo, const int* hi,
	          const amrex::Real* u, const int* u_lo, const int* u_hi,
		  const amrex::Real* q, const int* q_lo, const int* q_hi,
		  const amrex::Real* qaux, const int* qa_lo, const int* qa_hi, const int* idx);

  void ca_srctoprim(const int* lo,   const int* hi,
		    const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		    const amrex::Real* src, const int* src_lo, const int* src_hi,
		    const amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi, const int* idx);

  void ca_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_minimum_density
    (const amrex::Real* S_old, const int* s_old_lo, const int* s_old_hi,
     const amrex::Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const amrex::Real* vol, const int* vol_lo, const int* vol_hi,
     const int* lo, const int* hi,
     const amrex::Real* frac_change, const int* verbose, const int* idx);

  void ca_normalize_species
    (BL_FORT_FAB_ARG_3D(S_new), const int* lo, const int* hi, const int* idx);

  void ca_get_center(amrex::Real* center);
  void ca_set_center(amrex::Real* center);
  void ca_find_center(amrex::Real* data, amrex::Real* center, int* icen,
		      const amrex::Real* dx, const amrex::Real* problo);

  void ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new),
     const int& verbose, const int* idx);


  void ca_hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_xmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_ymomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_zmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_zerofill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_default_fill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  inline void ca_nullfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[])
#ifdef __cplusplus
  { amrex::Abort("Not supposed to call ca_nullfill"); }
#else
    ;
#endif

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

  void ca_check_timestep
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);

  void ca_ctu_update
    (const int* is_finest_level,
     const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     BL_FORT_FAB_ARG(q),
     BL_FORT_FAB_ARG(qaux),
     BL_FORT_FAB_ARG(src_q),
     BL_FORT_FAB_ARG(update),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG(pradial),
#endif
     D_DECL(BL_FORT_FAB_ARG(xarea),
            BL_FORT_FAB_ARG(yarea),
            BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     amrex::Real* cflLoc,
     const int&  verbose,
     amrex::Real& mass_lost,
     amrex::Real& xmom_lost,
     amrex::Real& ymom_lost,
     amrex::Real& zmom_lost,
     amrex::Real& eden_lost,
     amrex::Real& xang_lost,
     amrex::Real& yang_lost,
     amrex::Real& zang_lost);


  void ca_mol_single_stage
    (const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     BL_FORT_FAB_ARG(q),
     BL_FORT_FAB_ARG(qaux),
     BL_FORT_FAB_ARG(update),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG(pradial),
#endif
     D_DECL(BL_FORT_FAB_ARG(xarea),
            BL_FORT_FAB_ARG(yarea),
            BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     amrex::Real* cflLoc,
     const int&  verbose);

  void flush_output();

  void ca_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumsquared
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumlocmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);

  void ca_sumlocmass2d
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir1, const int& idir2);
  
  void ca_sumlocsquaredmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);
  
  void ca_sumproduct
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(f1), BL_FORT_FAB_ARG_3D(f2),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* s);

  void ca_tempdiffextrap
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(TempDiffTerm));

  void make_enthalpy
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
           BL_FORT_FAB_ARG_3D(enth));

  void ca_fill_temp_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_enth_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));
  
  void ca_fill_spec_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_first_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_secnd_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_compute_div_tau_u
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(div_tau_u),
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const int* coord_type);


#ifdef __cplusplus
}
#endif


BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);

BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(old_state),
     const BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);

#endif
