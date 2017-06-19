#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, int istage, int nstages)
{

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[istage];

  BL_PROFILE_VAR("Castro::construct_mol_hydro_source()", CA_HYDRO);

  MultiFab flux_mf[BL_SPACEDIM];
  MultiFab qe[BL_SPACEDIM];

  for (int i = 0; i < BL_SPACEDIM; ++i) {
      flux_mf[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
  }

  MultiFab q_mf;
  q_mf.define(grids, dmap, NQ, NUM_GROW);

  MultiFab flatn_mf;
  flatn_mf.define(grids, dmap, NQ, NUM_GROW);

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
	const Box& qbx = amrex::grow(bx, NUM_GROW);

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

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = k_stage[mfi];

	FArrayBox& vol = volume[mfi];

	q.resize(qbx, QVAR);
	qaux.resize(qbx, NQAUX);

	// convert the conservative state to the primitive variable state.
	// this fills both q and qaux.

	const int idx = mfi.tileIndex();

	Device::set_stream_index(idx);

	ca_ctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		   statein.dataPtr(), ARLIM_3D(statein.loVect()), ARLIM_3D(statein.hiVect()),
		   q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		   qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()));

	ca_mol_single_stage
	  (time,
	   lo, hi, domain_lo, domain_hi,
	   BL_TO_FORTRAN_3D(statein), 
	   BL_TO_FORTRAN_3D(stateout),
	   BL_TO_FORTRAN_3D(q),
	   BL_TO_FORTRAN_3D(flatn),
	   BL_TO_FORTRAN_3D(div),
	   BL_TO_FORTRAN_3D(qaux),
	   BL_TO_FORTRAN_3D(source_out),
	   dx, dt,
	   BL_TO_FORTRAN_3D(qe[0][mfi]),
	   BL_TO_FORTRAN_3D(qe[1][mfi]),
	   BL_TO_FORTRAN_3D(qe[2][mfi]),
	   BL_TO_FORTRAN_3D(qm),
	   BL_TO_FORTRAN_3D(qp),
	   BL_TO_FORTRAN_3D(sxm),
	   BL_TO_FORTRAN_3D(sxp),
	   BL_TO_FORTRAN_3D(sym),
	   BL_TO_FORTRAN_3D(syp),
	   BL_TO_FORTRAN_3D(szm),
	   BL_TO_FORTRAN_3D(szp),
	   BL_TO_FORTRAN_3D(flux_mf[0][mfi]),
	   BL_TO_FORTRAN_3D(flux_mf[1][mfi]),
	   BL_TO_FORTRAN_3D(flux_mf[2][mfi]),
	   BL_TO_FORTRAN_3D(area[0][mfi]),
	   BL_TO_FORTRAN_3D(area[1][mfi]),
	   BL_TO_FORTRAN_3D(area[2][mfi]),
	   BL_TO_FORTRAN_3D(volume[mfi]),
	   verbose);

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
