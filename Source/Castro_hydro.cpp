#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::cons_to_prim()
{

    BL_PROFILE("Castro::cons_to_prim()");

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(4);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]));

    }

}

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

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();

      // compute the flattening coefficient

      Array4<Real> const flatn_arr = flatn.array();

      ca_uflatten(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
                  BL_TO_FORTRAN_ANYD(q[mfi]),
                  BL_TO_FORTRAN_ANYD(flatn));

      for (int i = 0; i < 3; ++i) {
          ebx[i] = amrex::surroundingNodes(bx, i);
          gebx[i] = amrex::grow(ebx[i], 1);
      }

      tbx[0][0] = amrex::grow(bx, 1);
      tbx[0][1] = amrex::grow(ebx[0], 2, 1);
      tbx[0][2] = amrex::grow(ebx[0], 1, 1);
      tbx[1][0] = amrex::grow(ebx[1], 2, 1);
      tbx[1][1] = amrex::grow(bx, 1);
      tbx[1][2] = amrex::grow(ebx[1], 0, 1);
      tbx[2][0] = amrex::grow(ebx[2], 1, 1);
      tbx[2][1] = amrex::grow(ebx[2], 0, 1);
      tbx[2][2] = amrex::grow(bx, 1);

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

      ctu_ppm_states(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
                     AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                     BL_TO_FORTRAN_ANYD(q[mfi]),
                     BL_TO_FORTRAN_ANYD(flatn),
                     BL_TO_FORTRAN_ANYD(qaux[mfi]),
                     BL_TO_FORTRAN_ANYD(qm[0][0]),
                     BL_TO_FORTRAN_ANYD(qp[0][0]),
                     BL_TO_FORTRAN_ANYD(qm[1][1]),
                     BL_TO_FORTRAN_ANYD(qp[1][1]),
                     BL_TO_FORTRAN_ANYD(qm[2][2]),
                     BL_TO_FORTRAN_ANYD(qp[2][2]),
                     AMREX_ZFILL(dx), dt,
                     AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      // compute divu -- we'll use this later when doing the artifical viscosity
      divu(AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_ZFILL(dx),
           BL_TO_FORTRAN_ANYD(div));

      q_int.resize(obx, QVAR);
      Elixir elix_q_int = q_int.elixir();

      flux[0].resize(gebx[0], NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();

      qe[0].resize(gebx[0], NGDNV);
      Elixir elix_qe_x = qe[0].elixir();

      flux[1].resize(gebx[1], NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();

      qe[1].resize(gebx[1], NGDNV);
      Elixir elix_qe_y = qe[1].elixir();

      flux[2].resize(gebx[2], NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();

      qe[2].resize(gebx[2], NGDNV);
      Elixir elix_qe_z = qe[2].elixir();

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

      const amrex::Real hdt = 0.5*dt;

      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];
      const amrex::Real hdtdz = 0.5*dt/dx[2];

      const amrex::Real cdtdx = dt/dx[0]/3.0;
      const amrex::Real cdtdy = dt/dx[1]/3.0;
      const amrex::Real cdtdz = dt/dx[2]/3.0;

      // compute F^x
      // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cxbx = amrex::grow(ebx[0], IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cxbx.loVect()), AMREX_ARLIM_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[0][0]),
                          BL_TO_FORTRAN_ANYD(qp[0][0]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          1, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
      transx_on_ystates(AMREX_ARLIM_ANYD(tbx[1][0].loVect()), AMREX_ARLIM_ANYD(tbx[1][0].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[1][1]),
                        BL_TO_FORTRAN_ANYD(qm[1][0]),
                        BL_TO_FORTRAN_ANYD(qp[1][1]),
                        BL_TO_FORTRAN_ANYD(qp[1][0]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      transx_on_zstates(AMREX_ARLIM_ANYD(tbx[2][0].loVect()), AMREX_ARLIM_ANYD(tbx[2][0].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[2][2]),
                        BL_TO_FORTRAN_ANYD(qm[2][0]),
                        BL_TO_FORTRAN_ANYD(qp[2][2]),
                        BL_TO_FORTRAN_ANYD(qp[2][0]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ebx[1], IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cybx.loVect()), AMREX_ARLIM_ANYD(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[1][1]),
                          BL_TO_FORTRAN_ANYD(qp[1][1]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          2, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      transy_on_xstates(AMREX_ARLIM_ANYD(tbx[0][1].loVect()), AMREX_ARLIM_ANYD(tbx[0][1].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[0][0]),
                        BL_TO_FORTRAN_ANYD(qm[0][1]),
                        BL_TO_FORTRAN_ANYD(qp[0][0]),
                        BL_TO_FORTRAN_ANYD(qp[0][1]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      transy_on_zstates(AMREX_ARLIM_ANYD(tbx[2][1].loVect()), AMREX_ARLIM_ANYD(tbx[2][1].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[2][2]),
                        BL_TO_FORTRAN_ANYD(qm[2][1]),
                        BL_TO_FORTRAN_ANYD(qp[2][2]),
                        BL_TO_FORTRAN_ANYD(qp[2][1]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(ebx[2], IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(czbx.loVect()), AMREX_ARLIM_ANYD(czbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[2][2]),
                          BL_TO_FORTRAN_ANYD(qp[2][2]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          3, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      transz_on_xstates(AMREX_ARLIM_ANYD(tbx[0][2].loVect()), AMREX_ARLIM_ANYD(tbx[0][2].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[0][0]),
                        BL_TO_FORTRAN_ANYD(qm[0][2]),
                        BL_TO_FORTRAN_ANYD(qp[0][0]),
                        BL_TO_FORTRAN_ANYD(qp[0][2]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      transz_on_ystates(AMREX_ARLIM_ANYD(tbx[1][2].loVect()), AMREX_ARLIM_ANYD(tbx[1][2].hiVect()),
                        BL_TO_FORTRAN_ANYD(qm[1][1]),
                        BL_TO_FORTRAN_ANYD(qm[1][2]),
                        BL_TO_FORTRAN_ANYD(qp[1][1]),
                        BL_TO_FORTRAN_ANYD(qp[1][2]),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

      //
      // Use qx?, q?yz, q?zy to compute final x-flux
      //

      // compute F^{y|z}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cyzbx = amrex::grow(ebx[1], IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fyz
      // rftmp1 = rfyz
      // qgdnvtmp1 = qgdnvyz
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cyzbx.loVect()), AMREX_ARLIM_ANYD(cyzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[1][2]),
                          BL_TO_FORTRAN_ANYD(qp[1][2]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          2, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(ebx[2], IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(czybx.loVect()), AMREX_ARLIM_ANYD(czybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[2][1]),
                          BL_TO_FORTRAN_ANYD(qp[2][1]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          3, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

      transyz(AMREX_ARLIM_ANYD(ebx[0].loVect()), AMREX_ARLIM_ANYD(ebx[0].hiVect()),
              BL_TO_FORTRAN_ANYD(qm[0][0]),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qp[0][0]),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
              BL_TO_FORTRAN_ANYD(ftmp2),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdy, hdtdz);

      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(ebx[0].loVect()), AMREX_ARLIM_ANYD(ebx[0].hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          1, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(ebx[2], IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(czxbx.loVect()), AMREX_ARLIM_ANYD(czxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[2][0]),
                          BL_TO_FORTRAN_ANYD(qp[2][0]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          3, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(ebx[0], IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cxzbx.loVect()), AMREX_ARLIM_ANYD(cxzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[0][2]),
                          BL_TO_FORTRAN_ANYD(qp[0][2]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          1, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

      transxz(AMREX_ARLIM_ANYD(ebx[1].loVect()), AMREX_ARLIM_ANYD(ebx[1].hiVect()),
              BL_TO_FORTRAN_ANYD(qm[1][1]),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qp[1][1]),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp2),
              BL_TO_FORTRAN_ANYD(ftmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              hdt, hdtdx, hdtdz);

      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(ebx[1].loVect()), AMREX_ARLIM_ANYD(ebx[1].hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1]),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          2, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = amrex::grow(ebx[0], IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cxybx.loVect()), AMREX_ARLIM_ANYD(cxybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[0][1]),
                          BL_TO_FORTRAN_ANYD(qp[0][1]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          1, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ebx[1], IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(cyxbx.loVect()), AMREX_ARLIM_ANYD(cyxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qm[1][0]),
                          BL_TO_FORTRAN_ANYD(qp[1][0]), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          2, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      transxy(AMREX_ARLIM_ANYD(ebx[2].loVect()), AMREX_ARLIM_ANYD(ebx[2].hiVect()),
              BL_TO_FORTRAN_ANYD(qm[2][2]),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qp[2][2]),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
              BL_TO_FORTRAN_ANYD(ftmp2),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdx, hdtdy);

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      cmpflx_plus_godunov(AMREX_ARLIM_ANYD(ebx[2].loVect()), AMREX_ARLIM_ANYD(ebx[2].hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[2]),
                          BL_TO_FORTRAN_ANYD(q_int),
                          BL_TO_FORTRAN_ANYD(qe[2]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          3, AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi));

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          int idir_f = idir + 1;

          // Apply artificial viscosity to the fluxes.

          apply_av(AMREX_ARLIM_ANYD(nbx.loVect()), AMREX_ARLIM_ANYD(nbx.hiVect()),
                   idir_f, AMREX_ZFILL(dx),
                   BL_TO_FORTRAN_ANYD(div),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]));

          // Ensure species fluxes are normalized properly.

          normalize_species_fluxes(AMREX_ARLIM_ANYD(nbx.loVect()), AMREX_ARLIM_ANYD(nbx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(flux[idir]));

      }



      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();

      // conservative update

      ctu_consup(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
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

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

        scale_flux(AMREX_ARLIM_ANYD(nbx.loVect()), AMREX_ARLIM_ANYD(nbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

        // Store the fluxes from this advance.

        store_flux(AMREX_ARLIM_ANYD(nbx.loVect()), AMREX_ARLIM_ANYD(nbx.hiVect()),
                   BL_TO_FORTRAN_ANYD((*fluxes[idir])[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]));

      } // idir loop

    } // MFIter loop

  } // OMP loop

}
