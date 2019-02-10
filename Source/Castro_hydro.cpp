#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, int istage, int nstages)
{

  BL_PROFILE_VAR("Castro::construct_mol_hydro_source()", CA_HYDRO);

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[istage];

  MultiFab flux[BL_SPACEDIM];
  MultiFab qe[BL_SPACEDIM];

  for (int i = 0; i < BL_SPACEDIM; ++i) {
      flux[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
  }

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      const Box& qbx = mfi.growntilebox(NUM_GROW);

      // Convert the conservative state to the primitive variable state.
      // This fills both q and qaux.

      AsyncFab q(qbx, NQ);
      AsyncFab qaux(qbx, NQAUX);

#pragma gpu
      ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q.hostFab()),
                 BL_TO_FORTRAN_ANYD(qaux.hostFab()));

      const Box& obx = mfi.growntilebox(1);
      const Box& tbx = mfi.growntilebox(2);

      AsyncFab div(obx, 1);
      
      // Compute divergence of velocity field.

#pragma gpu
      ca_divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
              AMREX_REAL_ANYD(dx),
              BL_TO_FORTRAN_ANYD(q.hostFab()),
              BL_TO_FORTRAN_ANYD(div.hostFab()));

      AsyncFab flatn(obx, 1);

      // Compute flattening coefficient for slope calculations.
#pragma gpu
      ca_uflaten
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q.hostFab()),
           BL_TO_FORTRAN_ANYD(flatn.hostFab()));

      AsyncFab qm(tbx, 3*NQ);
      AsyncFab qp(tbx, 3*NQ);
      
      // Do PPM reconstruction to the zone edges.
#pragma gpu
      ca_ppm_reconstruct
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q.hostFab()),
           BL_TO_FORTRAN_ANYD(flatn.hostFab()),
           BL_TO_FORTRAN_ANYD(qm.hostFab()),
           BL_TO_FORTRAN_ANYD(qp.hostFab()));

      q.clear();
      flatn.clear();

      for (int idir = 0; idir < BL_SPACEDIM; ++idir) {

          const Box& ebx = mfi.nodaltilebox(idir);

          int idir_f = idir + 1;

#pragma gpu
          ca_construct_flux
              (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
               AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
               AMREX_REAL_ANYD(dx), dt,
               idir_f,
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(div.hostFab()),
               BL_TO_FORTRAN_ANYD(qaux.hostFab()),
               BL_TO_FORTRAN_ANYD(qm.hostFab()),
               BL_TO_FORTRAN_ANYD(qp.hostFab()),
               BL_TO_FORTRAN_ANYD(qe[idir][mfi]),
               BL_TO_FORTRAN_ANYD(flux[idir][mfi]),
               BL_TO_FORTRAN_ANYD(area[idir][mfi]));

          auto const flux_fab = (flux[idir]).array(mfi);
          auto       fluxes_fab = (*fluxes[idir]).array(mfi);
          const int numcomp = NUM_STATE;
          const Real scale = b_mol[istage];

          AMREX_HOST_DEVICE_FOR_4D(ebx, numcomp, i, j, k, n,
          {
              fluxes_fab(i,j,k,n) += scale * flux_fab(i,j,k,n);
          });

      }

      div.clear();
      qaux.clear();
      qm.clear();
      qp.clear();

#pragma gpu
      ca_construct_hydro_update
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_REAL_ANYD(dx), dt,
           b_mol[istage],
           BL_TO_FORTRAN_ANYD(qe[0][mfi]),
           BL_TO_FORTRAN_ANYD(qe[1][mfi]),
           BL_TO_FORTRAN_ANYD(qe[2][mfi]),
           BL_TO_FORTRAN_ANYD(flux[0][mfi]),
           BL_TO_FORTRAN_ANYD(flux[1][mfi]),
           BL_TO_FORTRAN_ANYD(flux[2][mfi]),
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
           BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(hydro_source[mfi]));

  } // MFIter loop

  BL_PROFILE_VAR_STOP(CA_HYDRO);

  // Flush Fortran output

  if (verbose)
    flush_output();

}
