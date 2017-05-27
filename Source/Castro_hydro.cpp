#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, int istage, int nstages)
{

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering hydro advance" << std::endl << std::endl;

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();
  Real courno    = -1.0e+200;

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[istage];

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    FArrayBox flux[BL_SPACEDIM];
#if (BL_SPACEDIM <= 2)
    FArrayBox pradial(Box::TheUnitBox(),1);
#endif
    FArrayBox q, qaux;

    int priv_nstep_fsp = -1;

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

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  pradial.resize(amrex::surroundingNodes(bx,0),1);
	}
#endif
	  
	ca_mol_single_stage
	  (&time,
	   lo, hi, domain_lo, domain_hi,
	   BL_TO_FORTRAN(statein), 
	   BL_TO_FORTRAN(stateout),
	   BL_TO_FORTRAN(q),
	   BL_TO_FORTRAN(qaux),
	   BL_TO_FORTRAN(source_out),
	   dx, &dt,
	   D_DECL(BL_TO_FORTRAN(flux[0]),
		  BL_TO_FORTRAN(flux[1]),
		  BL_TO_FORTRAN(flux[2])),
#if (BL_SPACEDIM < 3)
	   BL_TO_FORTRAN(pradial),
#endif
	   D_DECL(BL_TO_FORTRAN(area[0][mfi]),
		  BL_TO_FORTRAN(area[1][mfi]),
		  BL_TO_FORTRAN(area[2][mfi])),
#if (BL_SPACEDIM < 3)
	   BL_TO_FORTRAN(dLogArea[0][mfi]),
#endif
	   BL_TO_FORTRAN(volume[mfi]),
	   &cflLoc, verbose);

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < BL_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[istage], flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
	}

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].plus(pradial,mfi.nodaltilebox(0),0,0,1);
	}
#endif
      } // MFIter loop

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
    {
      courno = std::max(courno,cflLoc);
    }
  }  // end of omp parallel region

  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();

  if (courno > 1.0) {
    std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
    if (hard_cfl_limit == 1)
      amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
  }

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << std::endl << "... Leaving hydro advance" << std::endl << std::endl;

}
