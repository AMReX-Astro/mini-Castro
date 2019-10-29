#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, int istage, int nstages)
{
    BL_PROFILE("Castro::construct_mol_hydro_source()");

    // This constructs the hydrodynamic source (essentially the flux
    // divergence) using method of lines integration.  The output, as a
    // update to the state, is stored in the k_mol array of MultiFabs.

    int finest_level = parent->finestLevel();

    auto dx = geom.CellSizeArray();

    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& k_stage = *k_mol[istage];

    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();

    std::vector<amrex::Real> b_mol{0.5, 0.5};

    const Real update_scale_factor = b_mol[istage];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox q;
        FArrayBox qaux;
        FArrayBox div;
        FArrayBox flatn;
        FArrayBox qm;
        FArrayBox qp;
        FArrayBox flux[3];
        FArrayBox qe[3];

        for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi) {

            const Box& box = mfi.tilebox();

            const Box& qbx = amrex::grow(box, 4);

            FArrayBox& state_old = Sborder[mfi];
            auto state_old_arr = Sborder[mfi].array();

            // Convert the conservative state to the primitive variable state.
            // This fills both q and qaux.

            q.resize(qbx, QVAR);
            Elixir elix_q = q.elixir();
            auto q_arr = q.array();

            qaux.resize(qbx, NQAUX);
            Elixir elix_qaux = qaux.elixir();
            auto qaux_arr = qaux.array();

            CASTRO_LAUNCH_LAMBDA(qbx, lbx,
            {
                ca_ctoprim(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                           AMREX_ARR4_TO_FORTRAN_ANYD(state_old_arr),
                           AMREX_ARR4_TO_FORTRAN_ANYD(q_arr),
                           AMREX_ARR4_TO_FORTRAN_ANYD(qaux_arr));
            });

            const Box& obx = amrex::grow(box, 1);
            const Box& tbx = amrex::grow(box, 2);

            div.resize(obx, 1);
            Elixir elix_div = div.elixir();
            auto div_arr = div.array();

            // Compute divergence of velocity field.

            CASTRO_LAUNCH_LAMBDA(obx, lbx,
            {
                ca_divu(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                        AMREX_ZFILL(dx.data()),
                        AMREX_ARR4_TO_FORTRAN_ANYD(q_arr),
                        AMREX_ARR4_TO_FORTRAN_ANYD(div_arr));
            });

            flatn.resize(obx, 1);
            Elixir elix_flatn = flatn.elixir();
            auto flatn_arr = flatn.array();

            // Compute flattening coefficient for slope calculations.

            CASTRO_LAUNCH_LAMBDA(obx, lbx,
            {
                ca_uflaten
                    (AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     AMREX_ARR4_TO_FORTRAN_ANYD(q_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(flatn_arr));
            });

            qm.resize(tbx, 3*QVAR);
            Elixir elix_qm = qm.elixir();
            auto qm_arr = qm.array();

            qp.resize(tbx, 3*QVAR);
            Elixir elix_qp = qp.elixir();
            auto qp_arr = qp.array();

            // Do PPM reconstruction to the zone edges.

            CASTRO_LAUNCH_LAMBDA(obx, lbx,
            {
                ca_ppm_reconstruct
                    (AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     AMREX_ARR4_TO_FORTRAN_ANYD(q_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(flatn_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qm_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qp_arr));
            });

            q.clear();
            flatn.clear();

            flux[0].resize(amrex::surroundingNodes(box, 0), NUM_STATE);
            Elixir elix_flux_x = flux[0].elixir();
            auto flux_x_arr = flux[0].array();

            flux[1].resize(amrex::surroundingNodes(box, 1), NUM_STATE);
            Elixir elix_flux_y = flux[1].elixir();
            auto flux_y_arr = flux[1].array();

            flux[2].resize(amrex::surroundingNodes(box, 2), NUM_STATE);
            Elixir elix_flux_z = flux[2].elixir();
            auto flux_z_arr = flux[2].array();

            qe[0].resize(amrex::surroundingNodes(box, 0), NGDNV);
            Elixir elix_qe_x = qe[0].elixir();
            auto qe_x_arr = qe[0].array();

            qe[1].resize(amrex::surroundingNodes(box, 1), NGDNV);
            Elixir elix_qe_y = qe[1].elixir();
            auto qe_y_arr = qe[1].array();

            qe[2].resize(amrex::surroundingNodes(box, 2), NGDNV);
            Elixir elix_qe_z = qe[2].elixir();
            auto qe_z_arr = qe[2].array();

            for (int idir = 0; idir < 3; ++idir) {

                const Box& ebx = amrex::surroundingNodes(box, idir);

                int idir_f = idir + 1;

                auto flux_arr = flux[idir].array();
                auto fluxes_arr = fluxes[idir]->array(mfi);
                auto qe_arr   = qe[idir].array();
                auto area_arr = area[idir][mfi].array();

                CASTRO_LAUNCH_LAMBDA(ebx, lbx,
                {
                    ca_construct_flux
                        (AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                         AMREX_ARLIM_ANYD(domain_lo), AMREX_ARLIM_ANYD(domain_hi),
                         AMREX_ZFILL(dx.data()), dt,
                         idir_f,
                         update_scale_factor,
                         AMREX_ARR4_TO_FORTRAN_ANYD(state_old_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(div_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(qaux_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(qm_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(qp_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(qe_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(flux_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(fluxes_arr),
                         AMREX_ARR4_TO_FORTRAN_ANYD(area_arr));
                });

            }

            div.clear();
            qaux.clear();
            qm.clear();
            qp.clear();

            auto area_x_arr = area[0][mfi].array();
            auto area_y_arr = area[1][mfi].array();
            auto area_z_arr = area[2][mfi].array();

            auto volume_arr = volume[mfi].array();
            auto hydro_source_arr = hydro_source[mfi].array();

            CASTRO_LAUNCH_LAMBDA(box, lbx,
            {
                ca_construct_hydro_update
                    (AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     AMREX_ZFILL(dx.data()), dt,
                     update_scale_factor,
                     AMREX_ARR4_TO_FORTRAN_ANYD(qe_x_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qe_y_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(qe_z_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(flux_x_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(flux_y_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(flux_z_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(area_x_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(area_y_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(area_z_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(volume_arr),
                     AMREX_ARR4_TO_FORTRAN_ANYD(hydro_source_arr));
            });

        } // MFIter loop

    } // OpenMP loop

}
