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
  Real courno    = -1.0e+200;

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[istage];

  BL_PROFILE_VAR("Castro::construct_mol_hydro_source()", CA_HYDRO);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    FArrayBox flux[BL_SPACEDIM];

    FArrayBox q, qaux;

    Real cflLoc = -1.0e+200;

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();
	const Box& qbx = amrex::grow(bx, NUM_GROW);

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

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

	ca_ctoprim(ARLIM_3D(qbx.loVect()), ARLIM_3D(qbx.hiVect()),
		   statein.dataPtr(), ARLIM_3D(statein.loVect()), ARLIM_3D(statein.hiVect()),
		   q.dataPtr(), ARLIM_3D(q.loVect()), ARLIM_3D(q.hiVect()),
		   qaux.dataPtr(), ARLIM_3D(qaux.loVect()), ARLIM_3D(qaux.hiVect()), &idx);


	// Allocate fabs for fluxes
	for (int i = 0; i < BL_SPACEDIM ; i++)  {
	  const Box& bxtmp = amrex::surroundingNodes(bx,i);
	  flux[i].resize(bxtmp,NUM_STATE);
	}

	ca_mol_single_stage
	  (&time,
	   lo, hi, domain_lo, domain_hi,
	   BL_TO_FORTRAN_3D(statein), 
	   BL_TO_FORTRAN_3D(stateout),
	   BL_TO_FORTRAN_3D(q),
	   BL_TO_FORTRAN_3D(qaux),
	   BL_TO_FORTRAN_3D(source_out),
	   dx, &dt,
	   D_DECL(BL_TO_FORTRAN_3D(flux[0]),
		  BL_TO_FORTRAN_3D(flux[1]),
		  BL_TO_FORTRAN_3D(flux[2])),
	   D_DECL(BL_TO_FORTRAN_3D(area[0][mfi]),
		  BL_TO_FORTRAN_3D(area[1][mfi]),
		  BL_TO_FORTRAN_3D(area[2][mfi])),
	   BL_TO_FORTRAN_3D(volume[mfi]),
	   &cflLoc, verbose, &idx);

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < BL_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[istage], flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
	}

      } // MFIter loop

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
    {
      courno = std::max(courno,cflLoc);
    }
  }  // end of omp parallel region

  BL_PROFILE_VAR_STOP(CA_HYDRO);

  // Flush Fortran output

  if (verbose)
    flush_output();

  if (courno > 1.0)
    amrex::Abort("CFL is too high at this level");

}
