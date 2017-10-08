#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, int istage, int nstages)
{

  BL_PROFILE_VAR("Castro::construct_mol_hydro_source()", CA_HYDRO);

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSizeF();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[istage];

  MultiFab flux_mf[BL_SPACEDIM];
  MultiFab qe[BL_SPACEDIM];

  for (int i = 0; i < BL_SPACEDIM; ++i) {
      flux_mf[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
  }

  MultiFab q_mf;
  q_mf.define(grids, dmap, NQ, NUM_GROW);

  MultiFab flatn_mf;
  flatn_mf.define(grids, dmap, 1, NUM_GROW);

  MultiFab div_mf;
  div_mf.define(grids, dmap, 1, 1);

  MultiFab qm_mf;
  qm_mf.define(grids, dmap, 3*NQ, 1);

  MultiFab qp_mf;
  qp_mf.define(grids, dmap, 3*NQ, 1);

  MultiFab qaux_mf;
  qaux_mf.define(grids, dmap, NQAUX, NUM_GROW);

  MultiFab sxm_mf;
  sxm_mf.define(grids, dmap, NQ, 2);

  MultiFab sxp_mf;
  sxp_mf.define(grids, dmap, NQ, 2);

  MultiFab sym_mf;
  sym_mf.define(grids, dmap, NQ, 2);

  MultiFab syp_mf;
  syp_mf.define(grids, dmap, NQ, 2);

  MultiFab szm_mf;
  szm_mf.define(grids, dmap, NQ, 2);

  MultiFab szp_mf;
  szp_mf.define(grids, dmap, NQ, 2);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();
	const Box& qbx = mfi.registerBox(amrex::grow(bx, NUM_GROW));

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	FArrayBox &q        = q_mf[mfi];
	FArrayBox &flatn    = flatn_mf[mfi];
	FArrayBox &div      = div_mf[mfi];
	FArrayBox &qaux     = qaux_mf[mfi];
	FArrayBox &qm       = qm_mf[mfi];
	FArrayBox &qp       = qp_mf[mfi];
	FArrayBox &sxm      = sxm_mf[mfi];
	FArrayBox &sxp      = sxp_mf[mfi];
	FArrayBox &sym      = sym_mf[mfi];
	FArrayBox &syp      = syp_mf[mfi];
	FArrayBox &szm      = szm_mf[mfi];
	FArrayBox &szp      = szp_mf[mfi];
	FArrayBox &statein  = Sborder[mfi];
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_out = hydro_source[mfi];

	FArrayBox& vol = volume[mfi];

	// convert the conservative state to the primitive variable state.
	// this fills both q and qaux.

#ifdef AMREX_USE_DEVICE
        Device::prepare_for_launch(qbx.loVect(), qbx.hiVect());
#endif

	ca_ctoprim
          (BL_TO_FORTRAN_BOX(qbx),
	   BL_TO_FORTRAN_ANYD(statein),
           BL_TO_FORTRAN_ANYD(q),
           BL_TO_FORTRAN_ANYD(qaux));

        const Box& obx = mfi.registerBox(amrex::grow(bx, 1));

#ifdef AMREX_USE_DEVICE
        Device::prepare_for_launch(obx.loVect(), obx.hiVect());
#endif

        ca_prepare_for_fluxes
          (BL_TO_FORTRAN_BOX(obx),
	   dx, dt,
	   BL_TO_FORTRAN_ANYD(q),
	   BL_TO_FORTRAN_ANYD(qaux),
	   BL_TO_FORTRAN_ANYD(flatn),
	   BL_TO_FORTRAN_ANYD(div),
	   BL_TO_FORTRAN_ANYD(sxm),
	   BL_TO_FORTRAN_ANYD(sxp),
	   BL_TO_FORTRAN_ANYD(sym),
	   BL_TO_FORTRAN_ANYD(syp),
	   BL_TO_FORTRAN_ANYD(szm),
	   BL_TO_FORTRAN_ANYD(szp));

        ca_prepare_profile
          (BL_TO_FORTRAN_BOX(obx),
	   BL_TO_FORTRAN_ANYD(q),
	   BL_TO_FORTRAN_ANYD(qm),
	   BL_TO_FORTRAN_ANYD(qp),
	   BL_TO_FORTRAN_ANYD(sxm),
	   BL_TO_FORTRAN_ANYD(sxp),
	   BL_TO_FORTRAN_ANYD(sym),
	   BL_TO_FORTRAN_ANYD(syp),
	   BL_TO_FORTRAN_ANYD(szm),
	   BL_TO_FORTRAN_ANYD(szp));

        for (int idir = 0; idir < BL_SPACEDIM; ++idir) {

            int idir_f = idir + 1;

            const Box& ebx = mfi.registerBox(amrex::surroundingNodes(bx, idir));

#ifdef AMREX_USE_DEVICE
            Device::prepare_for_launch(ebx.loVect(), ebx.hiVect());
#endif

            ca_construct_flux
              (BL_TO_FORTRAN_BOX(ebx),
               domain_lo, domain_hi,
               dx, dt,
               idir_f,
               BL_TO_FORTRAN_ANYD(statein),
               BL_TO_FORTRAN_ANYD(div),
               BL_TO_FORTRAN_ANYD(qaux),
               BL_TO_FORTRAN_ANYD(qm),
               BL_TO_FORTRAN_ANYD(qp),
               BL_TO_FORTRAN_ANYD(qe[idir][mfi]),
               BL_TO_FORTRAN_ANYD(flux_mf[idir][mfi]),
               BL_TO_FORTRAN_ANYD(area[idir][mfi]));

        }

#ifdef AMREX_USE_DEVICE
        Device::prepare_for_launch(bx.loVect(), bx.hiVect());
#endif

	ca_construct_hydro_update
          (BL_TO_FORTRAN_BOX(bx),
           dx, dt,
	   b_mol[istage],
           BL_TO_FORTRAN_ANYD(qe[0][mfi]),
	   BL_TO_FORTRAN_ANYD(qe[1][mfi]),
	   BL_TO_FORTRAN_ANYD(qe[2][mfi]),
	   BL_TO_FORTRAN_ANYD(flux_mf[0][mfi]),
	   BL_TO_FORTRAN_ANYD(flux_mf[1][mfi]),
	   BL_TO_FORTRAN_ANYD(flux_mf[2][mfi]),
	   BL_TO_FORTRAN_ANYD(area[0][mfi]),
	   BL_TO_FORTRAN_ANYD(area[1][mfi]),
	   BL_TO_FORTRAN_ANYD(area[2][mfi]),
	   BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(source_out));

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < BL_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[istage], flux_mf[i][mfi], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
	}

      } // MFIter loop

  }  // end of omp parallel region

  BL_PROFILE_VAR_STOP(CA_HYDRO);

  // Flush Fortran output

  if (verbose)
    flush_output();

}
