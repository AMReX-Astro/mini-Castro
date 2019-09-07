#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void network_init();

  void network_finalize();

  void eos_init();

  void eos_finalize();

  void probinit_finalize();

  void ca_get_num_spec(int* nspec);

  void ca_get_spec_names(int* spec_names, int* ispec, int* len);

  void ca_get_qvar(int* qvar);

  AMREX_GPU_DEVICE
  void ca_ctoprim(const int* lo, const int* hi,
                  const amrex::Real* u, const int* u_lo, const int* u_hi,
                  const amrex::Real* q, const int* q_lo, const int* q_hi,
                  const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  void ca_enforce_minimum_density
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_normalize_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_compute_temp
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

  AMREX_GPU_DEVICE
  void ca_divu
    (const int* lo, const int* hi,
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(div));

  AMREX_GPU_DEVICE
  void ca_uflaten
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn));

  AMREX_GPU_DEVICE
  void ca_ppm_reconstruct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp));

  AMREX_GPU_DEVICE
  void ca_construct_flux
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real dt,
     const int idir,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     BL_FORT_FAB_ARG_3D(qe),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area));

  AMREX_GPU_DEVICE
  void ca_construct_hydro_update
    (const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real dt,
     const amrex::Real stage_weight,
     BL_FORT_FAB_ARG_3D(qx),
     BL_FORT_FAB_ARG_3D(qy),
     BL_FORT_FAB_ARG_3D(qz),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     BL_FORT_FAB_ARG_3D(xarea),
     BL_FORT_FAB_ARG_3D(yarea),
     BL_FORT_FAB_ARG_3D(zarea),
     const BL_FORT_FAB_ARG_3D(volume),
     BL_FORT_FAB_ARG_3D(update));

  void ca_initdata
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo, const amrex::Real* probhi);

  void ca_denerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     const amrex::Real* den, const int* den_lo, const int* den_hi,
     const int8_t tagval, const int8_t clearval);

  void calculate_blast_radius
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo, const amrex::Real* probhi,
     amrex::Real* blast_mass, amrex::Real* blast_radius,
     const amrex::Real max_density);

#ifdef __cplusplus
}
#endif



#endif
