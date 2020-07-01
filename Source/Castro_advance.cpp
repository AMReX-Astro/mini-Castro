
#include <Castro.H>
#include <Castro_F.H>

#include <cmath>
#include <climits>

using std::string;
using namespace amrex;

Real
Castro::advance (Real time, Real dt, int  amr_iteration, int  amr_ncycle)

  // the main driver for a single level.  This will do either the SDC
  // algorithm or the Strang-split reactions algorithm.
  //
  // arguments:
  //    time          : the current simulation time
  //    dt            : the timestep to advance (e.g., go from time to
  //                    time + dt)
  //    amr_iteration : where we are in the current AMR subcycle.  Each
  //                    level will take a number of steps to reach the
  //                    final time of the coarser level below it.  This
  //                    counter starts at 1
  //    amr_ncycle    : the number of subcycles at this level

{
    BL_PROFILE("Castro::advance()");

    // Swap the new data from the last timestep into the old state data.

    state[0].allocOldData();
    state[0].swapTimeLevels(dt);

    // Allocate space for the MultiFabs we need during the step.

    hydro_source.define(grids,dmap,NUM_STATE,0);
    Sborder.define(grids, dmap, NUM_STATE, 4);

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Clean the old-time state, in case we came in after a regrid
    // where the state could be thermodynamically inconsistent.

    clean_state(S_old);

    // Initialize the new-time data from the old time data.

    for (MFIter mfi(S_old, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        auto S_old_arr = S_old[mfi].array();
        auto S_new_arr = S_new[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            copy(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                 AMREX_ARR4_TO_FORTRAN_ANYD(S_new_arr),
                 AMREX_ARR4_TO_FORTRAN_ANYD(S_old_arr));
        });
    }

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we first do a local copy to an array Sborder,
    // then we do a halo exchange to fill its ghost zones.

    for (MFIter mfi(S_old, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        auto S_old_arr = S_old[mfi].array();
        auto Sbord_arr = Sborder[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            copy(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                 AMREX_ARR4_TO_FORTRAN_ANYD(Sbord_arr),
                 AMREX_ARR4_TO_FORTRAN_ANYD(S_old_arr));
        });
    }

    Sborder.FillBoundary(geom.periodicity());

    // Make the temporarily expanded state thermodynamically consistent after the fill.

    clean_state(Sborder);

    // Construct the hydro source.

    construct_hydro_source(dt);

    // Add it to the state, scaled by the timestep.

    for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        auto S_new_arr = S_new[mfi].array();
        auto hydro_arr = hydro_source[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            axpy(AMREX_INT_ANYD(lbx.loVect()), AMREX_INT_ANYD(lbx.hiVect()),
                 AMREX_ARR4_TO_FORTRAN_ANYD(S_new_arr),
                 AMREX_ARR4_TO_FORTRAN_ANYD(hydro_arr),
                 dt);
        });
    }

    // Make the state thermodynamically consistent.

    clean_state(S_new);

    // Clear our temporary MultiFabs.

    hydro_source.clear();
    Sborder.clear();

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts() / getLevel(0).grids.numPts();

    return dt;
}
