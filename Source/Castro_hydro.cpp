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

  hydro_source.setVal(0.0);

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).

    FArrayBox flatn;
    FArrayBox dq;
    FArrayBox q;
    FArrayBox qaux;

    // The terms of qm and qp with i == j are the edge states that
    // come out of the PPM edge state prediction. The terms with
    // i /= j include transverse corrections.

    FArrayBox qm[3][3], qp[3][3];

    amrex::Array<Box, 3> ebx;
    amrex::Array<Box, 3> gebx;
    amrex::Array<amrex::Array<Box, 3>, 3> tbx;

    FArrayBox div;
    FArrayBox q_int;
    FArrayBox ftmp1, ftmp2;
    FArrayBox qgdnvtmp1, qgdnvtmp2;
    FArrayBox ql, qr;
    FArrayBox flux[3], qe[3];
    FArrayBox pdivu;

    for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      const Box& qbx = mfi.growntilebox(4);

      q.resize(qbx, QVAR);
      Elixir elix_q = q.elixir();

      qaux.resize(qbx, NQAUX);
      Elixir elix_qaux = qaux.elixir();

      // Convert the conservative state to the primitive variable state.

      ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q),
                 BL_TO_FORTRAN_ANYD(qaux));

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();

      // compute the flattening coefficient

      Array4<Real> const flatn_arr = flatn.array();

      ca_uflatten(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
                  BL_TO_FORTRAN_ANYD(q),
                  BL_TO_FORTRAN_ANYD(flatn));

      Elixir elix_flux[3];
      Elixir elix_qe[3];

      for (int i = 0; i < 3; ++i) {
          ebx[i] = amrex::surroundingNodes(bx, i);
          gebx[i] = amrex::grow(ebx[i], 1);

          flux[i].resize(gebx[i], NUM_STATE);
          elix_flux[i] = flux[i].elixir();

          qe[i].resize(gebx[i], NGDNV);
          elix_qe[i] = qe[i].elixir();
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

      Elixir elix_qm[3][3];
      Elixir elix_qp[3][3];

      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
              qm[i][j].resize(tbx[i][j], QVAR);
              elix_qm[i][j] = qm[i][j].elixir();

              qp[i][j].resize(obx, QVAR);
              elix_qp[i][j] = qp[i][j].elixir();
          }
      }

      int idir, idir_f;
      int idir_t1, idir_t1_f;
      int idir_t2, idir_t2_f;

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      // compute divu -- we'll use this later when doing the artificial viscosity
      divu(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q),
           AMREX_ZFILL(dx),
           BL_TO_FORTRAN_ANYD(div));

      q_int.resize(obx, QVAR);
      Elixir elix_q_int = q_int.elixir();

      ftmp1.resize(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1.elixir();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();

      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();

      ql.resize(obx, QVAR);
      Elixir elix_ql = ql.elixir();

      qr.resize(obx, QVAR);
      Elixir elix_qr = qr.elixir();

      const amrex::Real hdtdx[3] = {0.5*dt/dx[0], 0.5*dt/dx[1], 0.5*dt/dx[2]};
      const amrex::Real cdtdx[3] = {dt/dx[0]/3.0, dt/dx[1]/3.0, dt/dx[2]/3.0};

      for (idir = 0; idir < 3; ++idir) {

          idir_f = idir + 1;

          trace_ppm(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
                    AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                    idir_f,
                    BL_TO_FORTRAN_ANYD(q),
                    BL_TO_FORTRAN_ANYD(qaux),
                    BL_TO_FORTRAN_ANYD(flatn),
                    BL_TO_FORTRAN_ANYD(qm[idir][idir]),
                    BL_TO_FORTRAN_ANYD(qp[idir][idir]),
                    AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi),
                    AMREX_ZFILL(dx), dt);

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

          cmpflx_plus_godunov(AMREX_ARLIM_ANYD(tbx[idir][idir].loVect()), AMREX_ARLIM_ANYD(tbx[idir][idir].hiVect()),
                              BL_TO_FORTRAN_ANYD(qm[idir][idir]),
                              BL_TO_FORTRAN_ANYD(qp[idir][idir]), 1, 1,
                              BL_TO_FORTRAN_ANYD(ftmp1),
                              BL_TO_FORTRAN_ANYD(q_int),
                              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                              BL_TO_FORTRAN_ANYD(qaux),
                              idir_f, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

          trans1(AMREX_ARLIM_ANYD(tbx[idir_t1][idir].loVect()), AMREX_ARLIM_ANYD(tbx[idir_t1][idir].hiVect()),
                 idir_f, idir_t1_f,
                 BL_TO_FORTRAN_ANYD(qm[idir_t1][idir_t1]),
                 BL_TO_FORTRAN_ANYD(qm[idir_t1][idir]),
                 BL_TO_FORTRAN_ANYD(qp[idir_t1][idir_t1]),
                 BL_TO_FORTRAN_ANYD(qp[idir_t1][idir]),
                 BL_TO_FORTRAN_ANYD(qaux),
                 BL_TO_FORTRAN_ANYD(ftmp1),
                 BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                 cdtdx[idir]);

          trans1(AMREX_ARLIM_ANYD(tbx[idir_t2][idir].loVect()), AMREX_ARLIM_ANYD(tbx[idir_t2][idir].hiVect()),
                 idir_f, idir_t2_f,
                 BL_TO_FORTRAN_ANYD(qm[idir_t2][idir_t2]),
                 BL_TO_FORTRAN_ANYD(qm[idir_t2][idir]),
                 BL_TO_FORTRAN_ANYD(qp[idir_t2][idir_t2]),
                 BL_TO_FORTRAN_ANYD(qp[idir_t2][idir]),
                 BL_TO_FORTRAN_ANYD(qaux),
                 BL_TO_FORTRAN_ANYD(ftmp1),
                 BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                 cdtdx[idir]);

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

          // compute F^{1|2}
          cmpflx_plus_godunov(AMREX_ARLIM_ANYD(tbx[idir_t1][idir_t2].loVect()), AMREX_ARLIM_ANYD(tbx[idir_t1][idir_t2].hiVect()),
                              BL_TO_FORTRAN_ANYD(qm[idir_t1][idir_t2]),
                              BL_TO_FORTRAN_ANYD(qp[idir_t1][idir_t2]), 1, 1,
                              BL_TO_FORTRAN_ANYD(ftmp1),
                              BL_TO_FORTRAN_ANYD(q_int),
                              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                              BL_TO_FORTRAN_ANYD(qaux),
                              idir_t1_f, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

          // compute F^{2|1}
          cmpflx_plus_godunov(AMREX_ARLIM_ANYD(tbx[idir_t2][idir_t1].loVect()), AMREX_ARLIM_ANYD(tbx[idir_t2][idir_t1].hiVect()),
                              BL_TO_FORTRAN_ANYD(qm[idir_t2][idir_t1]),
                              BL_TO_FORTRAN_ANYD(qp[idir_t2][idir_t1]), 1, 1,
                              BL_TO_FORTRAN_ANYD(ftmp2),
                              BL_TO_FORTRAN_ANYD(q_int),
                              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                              BL_TO_FORTRAN_ANYD(qaux),
                              idir_t2_f, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

          // Compute the corrected idir interface states and fluxes
          trans2(AMREX_ARLIM_ANYD(ebx[idir].loVect()), AMREX_ARLIM_ANYD(ebx[idir].hiVect()),
                 idir_f, idir_t1_f, idir_t2_f,
                 BL_TO_FORTRAN_ANYD(qm[idir][idir]),
                 BL_TO_FORTRAN_ANYD(ql),
                 BL_TO_FORTRAN_ANYD(qp[idir][idir]),
                 BL_TO_FORTRAN_ANYD(qr),
                 BL_TO_FORTRAN_ANYD(qaux),
                 BL_TO_FORTRAN_ANYD(ftmp1),
                 BL_TO_FORTRAN_ANYD(ftmp2),
                 BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                 BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                 hdtdx[idir], hdtdx[idir_t1], hdtdx[idir_t2]);

          // Compute the final F^idir
          cmpflx_plus_godunov(AMREX_ARLIM_ANYD(ebx[idir].loVect()), AMREX_ARLIM_ANYD(ebx[idir].hiVect()),
                              BL_TO_FORTRAN_ANYD(ql),
                              BL_TO_FORTRAN_ANYD(qr), 1, 1,
                              BL_TO_FORTRAN_ANYD(flux[idir]),
                              BL_TO_FORTRAN_ANYD(q_int),
                              BL_TO_FORTRAN_ANYD(qe[idir]),
                              BL_TO_FORTRAN_ANYD(qaux),
                              idir_f, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      }

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          int idir_f = idir + 1;

          // Apply artificial viscosity to the fluxes.

          apply_av(AMREX_ARLIM_ANYD(ebx[idir].loVect()), AMREX_ARLIM_ANYD(ebx[idir].hiVect()),
                   idir_f, AMREX_ZFILL(dx),
                   BL_TO_FORTRAN_ANYD(div),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]));

          // Ensure species fluxes are normalized properly.

          normalize_species_fluxes(AMREX_ARLIM_ANYD(ebx[idir].loVect()), AMREX_ARLIM_ANYD(ebx[idir].hiVect()),
                                   BL_TO_FORTRAN_ANYD(flux[idir]));

          // Store the fluxes from this advance; we'll use these in
          // the flux register for doing the coarse-fine level sync.
          // The flux is scaled by dt * dA.

          store_flux(AMREX_ARLIM_ANYD(ebx[idir].loVect()), AMREX_ARLIM_ANYD(ebx[idir].hiVect()),
                     BL_TO_FORTRAN_ANYD((*fluxes[idir])[mfi]),
                     BL_TO_FORTRAN_ANYD(flux[idir]),
                     BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                     dt);

      }

      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();

      // conservative update

      ctu_consup(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q),
                 BL_TO_FORTRAN_ANYD(hydro_source[mfi]),
                 BL_TO_FORTRAN_ANYD(flux[0]),
                 BL_TO_FORTRAN_ANYD(flux[1]),
                 BL_TO_FORTRAN_ANYD(flux[2]),
                 BL_TO_FORTRAN_ANYD(qe[0]),
                 BL_TO_FORTRAN_ANYD(qe[1]),
                 BL_TO_FORTRAN_ANYD(qe[2]),
                 BL_TO_FORTRAN_ANYD(area[0][mfi]),
                 BL_TO_FORTRAN_ANYD(area[1][mfi]),
                 BL_TO_FORTRAN_ANYD(area[2][mfi]),
                 BL_TO_FORTRAN_ANYD(volume[mfi]),
                 BL_TO_FORTRAN_ANYD(pdivu),
                 AMREX_ZFILL(dx), dt);

    } // MFIter loop

  } // OMP loop

}
