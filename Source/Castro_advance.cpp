
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

    MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.

    AmrLevel::FillPatch(*this, Sborder, 4, time, State_Type, 0, NUM_STATE);

    // Make the temporarily expanded state thermodynamically consistent after the fill.

    clean_state(Sborder);

    // Construct the hydro source.

    construct_hydro_source(dt);

    // Add it to the state, scaled by the timestep.

    MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NUM_STATE, 0);

    // Make the state thermodynamically consistent.

    clean_state(S_new);

    // Clear our temporary MultiFabs.

    hydro_source.clear();
    Sborder.clear();

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts() / getLevel(0).grids.numPts();

    return dt;
}
