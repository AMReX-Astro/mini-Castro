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

  CASTRO_DEVICE
  void ca_ctoprim(const int* lo, const int* hi,
                  const amrex::Real* u, const int* u_lo, const int* u_hi,
                  const amrex::Real* q, const int* q_lo, const int* q_hi,
                  const amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  CASTRO_DEVICE
  void ca_enforce_minimum_density
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  CASTRO_DEVICE
  void ca_normalize_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  CASTRO_DEVICE
  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  CASTRO_DEVICE
  void ca_compute_temp
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  CASTRO_DEVICE
  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

  CASTRO_DEVICE
  void ca_divu
    (const int* lo, const int* hi,
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(div));

  CASTRO_DEVICE
  void ca_uflatten
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn));

  CASTRO_DEVICE
  void ca_ppm_reconstruct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp));

  CASTRO_DEVICE
  void ca_construct_flux
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real dt,
     const int idir,
     const amrex::Real update_scale_factor,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     BL_FORT_FAB_ARG_3D(qe),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(fluxes),
     BL_FORT_FAB_ARG_3D(area));

  CASTRO_DEVICE
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

  void scale_flux
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);
  void ctu_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(Sborder),
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(hydro_source),
     const BL_FORT_FAB_ARG_3D(flux0),
     const BL_FORT_FAB_ARG_3D(flux1),
     const BL_FORT_FAB_ARG_3D(flux2),
     const BL_FORT_FAB_ARG_3D(qe0),
     const BL_FORT_FAB_ARG_3D(qe1),
     const BL_FORT_FAB_ARG_3D(qe2),
     const BL_FORT_FAB_ARG_3D(area0),
     const BL_FORT_FAB_ARG_3D(area1),
     const BL_FORT_FAB_ARG_3D(area2),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(pdivu),
     const amrex::Real* dx, const amrex::Real dt);

  void normalize_species_fluxes
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(flux));

  void apply_av
    (const int* lo, const int* hi,
     const int idir, const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(Sborder),
     BL_FORT_FAB_ARG_3D(flux));

  void transx_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
                         const BL_FORT_FAB_ARG_3D(qx),
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transx_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
                         const BL_FORT_FAB_ARG_3D(qx),
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transy_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transy_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transz_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transz_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transyz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fyz),
               const BL_FORT_FAB_ARG_3D(fzy),
               const BL_FORT_FAB_ARG_3D(qy),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdy, const amrex::Real cdtdz);

  void transxz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxz),
               const BL_FORT_FAB_ARG_3D(fzx),
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdz);

  void transxy(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxy),
               const BL_FORT_FAB_ARG_3D(fyx),
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qy),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdy);

  void cmpflx_plus_godunov
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(ql),
     const BL_FORT_FAB_ARG_3D(qr), const int nc, const int comp,
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(q_int),
     BL_FORT_FAB_ARG_3D(qe),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(shk),
     const int idir, const int* domlo, const int* domhi);

  void divu
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div));

  void ctu_ppm_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      BL_FORT_FAB_ARG_3D(qxm),
      BL_FORT_FAB_ARG_3D(qxp),
      BL_FORT_FAB_ARG_3D(qym),
      BL_FORT_FAB_ARG_3D(qyp),
      BL_FORT_FAB_ARG_3D(qzm),
      BL_FORT_FAB_ARG_3D(qzp),
      const amrex::Real* dx, const amrex::Real dt,
      const int* domlo, const int* domhi);

  CASTRO_DEVICE
  void ca_initdata
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo, const amrex::Real* probhi);

  CASTRO_DEVICE
  void ca_denerror
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     const amrex::Real* den, const int* den_lo, const int* den_hi,
     const int8_t tagval, const int8_t clearval);

  CASTRO_DEVICE
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
