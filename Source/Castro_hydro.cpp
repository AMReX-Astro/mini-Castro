#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::construct_hydro_source(Real dt)
{

  BL_PROFILE("Castro::construct_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  int finest_level = parent->finestLevel();

  auto dx = geom.CellSizeArray();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      const Box& qbx = amrex::grow(bx, 4);

      Array4<Real> const state = Sborder[mfi].array();
      Array4<Real> const source = hydro_source[mfi].array();
      Array4<Real> const ar[3] = {area[0][mfi].array(), area[1][mfi].array(), area[2][mfi].array()};
      Array4<Real> const vol = volume[mfi].array();

      amrex::Array<Box, 3> ebx;
      amrex::Array<Box, 3> gebx;
      amrex::Array<amrex::Array<Box, 3>, 3> tbx;

      // Declare local storage now.
      // We apply an Elixir to ensure that memory is saved until it is no
      // longer needed (only relevant for the asynchronous case, usually on GPUs).

      FArrayBox q_fab(qbx, QVAR);
      Elixir elix_q = q_fab.elixir();
      Array4<Real> const q = q_fab.array();

      FArrayBox qaux_fab(qbx, NQAUX);
      Elixir elix_qaux = qaux_fab.elixir();
      Array4<Real> const qaux = qaux_fab.array();

      // Convert the conservative state to the primitive variable state.

      CASTRO_LAUNCH_LAMBDA(qbx, lbx,
      {
          ctoprim(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                  AMREX_ARR4_TO_FORTRAN_ANYD(state),
                  AMREX_ARR4_TO_FORTRAN_ANYD(q),
                  AMREX_ARR4_TO_FORTRAN_ANYD(qaux));
      });

      FArrayBox flux_fab[3], qe_fab[3];

      Elixir elix_flux[3];
      Elixir elix_qe[3];

      Array4<Real> flux[3];
      Array4<Real> qe[3];

      for (int i = 0; i < 3; ++i) {
          ebx[i] = amrex::surroundingNodes(bx, i);
          gebx[i] = amrex::grow(ebx[i], 1);

          flux_fab[i].resize(gebx[i], NUM_STATE);
          elix_flux[i] = flux_fab[i].elixir();
          flux[i] = flux_fab[i].array();

          qe_fab[i].resize(gebx[i], NGDNV);
          elix_qe[i] = qe_fab[i].elixir();
          qe[i] = qe_fab[i].array();
      }

      tbx[0][0] = amrex::grow(ebx[0], IntVect(0,1,1));
      tbx[0][1] = amrex::grow(ebx[0], IntVect(0,0,1));
      tbx[0][2] = amrex::grow(ebx[0], IntVect(0,1,0));
      tbx[1][0] = amrex::grow(ebx[1], IntVect(0,0,1));
      tbx[1][1] = amrex::grow(ebx[1], IntVect(1,0,1));
      tbx[1][2] = amrex::grow(ebx[1], IntVect(1,0,0));
      tbx[2][0] = amrex::grow(ebx[2], IntVect(0,1,0));
      tbx[2][1] = amrex::grow(ebx[2], IntVect(1,0,0));
      tbx[2][2] = amrex::grow(ebx[2], IntVect(1,1,0));
      
      // The terms of qm and qp with i == j are the edge states that
      // come out of the PPM edge state prediction. The terms with
      // i /= j include transverse corrections.

      FArrayBox qm_fab[3][3], qp_fab[3][3];

      Elixir elix_qm[3][3];
      Elixir elix_qp[3][3];

      Array4<Real> qm[3][3];
      Array4<Real> qp[3][3];

      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              qm_fab[i][j].resize(tbx[i][j], QVAR);
              elix_qm[i][j] = qm_fab[i][j].elixir();
              qm[i][j] = qm_fab[i][j].array();

              qp_fab[i][j].resize(tbx[i][j], QVAR);
              elix_qp[i][j] = qp_fab[i][j].elixir();
              qp[i][j] = qp_fab[i][j].array();
          }
      }

      int idir, idir_f;
      int idir_t1, idir_t1_f;
      int idir_t2, idir_t2_f;

      FArrayBox div_fab(obx, 1);
      Elixir elix_div = div_fab.elixir();
      Array4<Real> const div = div_fab.array();

      // Compute divu -- we'll use this later when doing the artificial viscosity
      CASTRO_LAUNCH_LAMBDA(obx, lbx,
      {
          divu(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
               AMREX_ARR4_TO_FORTRAN_ANYD(q),
               AMREX_ZFILL(dx.data()),
               AMREX_ARR4_TO_FORTRAN_ANYD(div));
      });

      FArrayBox q_int_fab(obx, QVAR);
      Elixir elix_q_int = q_int_fab.elixir();
      Array4<Real> const q_int = q_int_fab.array();

      FArrayBox ftmp1_fab(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1_fab.elixir();
      Array4<Real> const ftmp1 = ftmp1_fab.array();

      FArrayBox ftmp2_fab(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2_fab.elixir();
      Array4<Real> const ftmp2 = ftmp2_fab.array();

      FArrayBox qgdnvtmp1_fab(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1_fab.elixir();
      Array4<Real> const qgdnvtmp1 = qgdnvtmp1_fab.array();

      FArrayBox qgdnvtmp2_fab(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2_fab.elixir();
      Array4<Real> const qgdnvtmp2 = qgdnvtmp2_fab.array();

      FArrayBox ql_fab(obx, QVAR);
      Elixir elix_ql = ql_fab.elixir();
      Array4<Real> const ql = ql_fab.array();

      FArrayBox qr_fab(obx, QVAR);
      Elixir elix_qr = qr_fab.elixir();
      Array4<Real> const qr = qr_fab.array();

      const amrex::Real hdtdx[3] = {0.5*dt/dx[0], 0.5*dt/dx[1], 0.5*dt/dx[2]};
      const amrex::Real cdtdx[3] = {dt/dx[0]/3.0, dt/dx[1]/3.0, dt/dx[2]/3.0};

      for (idir = 0; idir < 3; ++idir) {

          idir_f = idir + 1;

          CASTRO_LAUNCH_LAMBDA(obx, lbx,
          {
              trace_ppm(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                        AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                        idir_f,
                        AMREX_ARR4_TO_FORTRAN_ANYD(q),
                        AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                        AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir][idir]),
                        AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir][idir]),
                        AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi),
                        AMREX_ZFILL(dx.data()), dt);
          });

      }

      for (idir = 0; idir < 3; ++idir) {

          if (idir == 0) {
              idir_t1 = 1;
              idir_t2 = 2;
          }
          else if (idir == 1) {
              idir_t1 = 0;
              idir_t2 = 2;
          }
          else {
              idir_t1 = 0;
              idir_t2 = 1;
          }

          idir_f = idir + 1;
          idir_t1_f = idir_t1 + 1;
          idir_t2_f = idir_t2 + 1;

          // Compute the flux in this coordinate direction
          CASTRO_LAUNCH_LAMBDA(tbx[idir][idir], lbx,
          {
              compute_flux(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir][idir]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir][idir]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(ftmp1),
                           AMREX_ARR4_TO_FORTRAN_ANYD(q_int),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp1),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                           idir_f);
          });

          // Update the states in one of the two orthogonal directions using the
          // transverse flux direction in this coordinate direction.
          CASTRO_LAUNCH_LAMBDA(tbx[idir_t1][idir], lbx,
          {
              trans1(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     idir_f, idir_t1_f,
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t1][idir_t1]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t1][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t1][idir_t1]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t1][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                     AMREX_ARR4_TO_FORTRAN_ANYD(ftmp1),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp1),
                     cdtdx[idir]);
          });

          // Do the same for the other orthogonal direction.
          CASTRO_LAUNCH_LAMBDA(tbx[idir_t2][idir], lbx,
          {
              trans1(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     idir_f, idir_t2_f,
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t2][idir_t2]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t2][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t2][idir_t2]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t2][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                     AMREX_ARR4_TO_FORTRAN_ANYD(ftmp1),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp1),
                     cdtdx[idir]);
          });

      }

      for (idir = 0; idir < 3; ++idir) {

          if (idir == 0) {
              idir_t1 = 1;
              idir_t2 = 2;
          }
          else if (idir == 1) {
              idir_t1 = 0;
              idir_t2 = 2;
          }
          else {
              idir_t1 = 0;
              idir_t2 = 1;
          }

          idir_f = idir + 1;
          idir_t1_f = idir_t1 + 1;
          idir_t2_f = idir_t2 + 1;

          // Compute F^{1|2}, the flux in direction 1 given the transverse flux correction
          // from direction 2.
          CASTRO_LAUNCH_LAMBDA(tbx[idir_t1][idir_t2], lbx,
          {
              compute_flux(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t1][idir_t2]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t1][idir_t2]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(ftmp1),
                           AMREX_ARR4_TO_FORTRAN_ANYD(q_int),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp1),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                           idir_t1_f);
          });

          // Compute F^{2|1}, the flux in direction 2 given the transverse flux correction
          // from direction 1.
          CASTRO_LAUNCH_LAMBDA(tbx[idir_t2][idir_t1], lbx,
          {                               
              compute_flux(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir_t2][idir_t1]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir_t2][idir_t1]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(ftmp2),
                           AMREX_ARR4_TO_FORTRAN_ANYD(q_int),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp2),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                           idir_t2_f);
          });

          // Compute the corrected idir interface states, given the two transverse fluxes.
          CASTRO_LAUNCH_LAMBDA(ebx[idir], lbx,
          {
              trans2(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     idir_f, idir_t1_f, idir_t2_f,
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm[idir][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(ql),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp[idir][idir]),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                     AMREX_ARR4_TO_FORTRAN_ANYD(ftmp1),
                     AMREX_ARR4_TO_FORTRAN_ANYD(ftmp2),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp1),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qgdnvtmp2),
                     hdtdx[idir], hdtdx[idir_t1], hdtdx[idir_t2]);
          });

          // Compute the final flux in direction idir, given the corrected interface states.
          CASTRO_LAUNCH_LAMBDA(ebx[idir], lbx,
          {
              compute_flux(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                           AMREX_ARR4_TO_FORTRAN_ANYD(ql),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qr),
                           AMREX_ARR4_TO_FORTRAN_ANYD(flux[idir]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(q_int),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qe[idir]),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qaux),
                           idir_f);
          });

      }

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          int idir_f = idir + 1;

          // Apply artificial viscosity to the fluxes.

          CASTRO_LAUNCH_LAMBDA(ebx[idir], lbx,
          {
              apply_av(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                       idir_f, AMREX_ZFILL(dx.data()),
                       AMREX_ARR4_TO_FORTRAN_ANYD(div),
                       AMREX_ARR4_TO_FORTRAN_ANYD(state),
                       AMREX_ARR4_TO_FORTRAN_ANYD(flux[idir]));
          });

          // Ensure species fluxes are normalized properly.

          CASTRO_LAUNCH_LAMBDA(ebx[idir], lbx,
          {
              normalize_species_fluxes(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                       AMREX_ARR4_TO_FORTRAN_ANYD(flux[idir]));
          });

      }

      // Construct the conservative update source term.

      CASTRO_LAUNCH_LAMBDA(bx, lbx,
      {
          fill_hydro_source(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                            AMREX_ARR4_TO_FORTRAN_ANYD(state),
                            AMREX_ARR4_TO_FORTRAN_ANYD(q),
                            AMREX_ARR4_TO_FORTRAN_ANYD(source),
                            AMREX_ARR4_TO_FORTRAN_ANYD(flux[0]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(flux[1]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(flux[2]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(qe[0]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(qe[1]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(qe[2]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(ar[0]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(ar[1]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(ar[2]),
                            AMREX_ARR4_TO_FORTRAN_ANYD(vol),
                            AMREX_ZFILL(dx.data()), dt);
      });

  } // MFIter loop

}
