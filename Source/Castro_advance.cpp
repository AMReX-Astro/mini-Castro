
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
    q.define(grids, dmap, QVAR, 4);
    qaux.define(grids, dmap, NQAUX, 4);

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	fluxes[dir]->setVal(0.0);

    // Do the advance.

    do_advance(time, dt);

    // Update the flux registers.

    FluxRegCrseInit();
    FluxRegFineAdd();

    // Clear our temporary arrays.

    hydro_source.clear();
    q.clear();
    qaux.clear();

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts() / getLevel(0).grids.numPts();

    return dt;
}

Real
Castro::do_advance (Real time, Real dt)
{

    // this routine will advance the old state data (called S_old here)
    // to the new time, for a single level.  The new data is called
    // S_new here.  The update includes reactions (if we are not doing
    // SDC), hydro, and the source terms.

    BL_PROFILE("Castro::do_advance()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Clean the old-time state, in case we came in after a regrid
    // where the state could be thermodynamically inconsistent.
    clean_state(S_old);

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.
    Sborder.define(grids, dmap, NUM_STATE, 4);
    AmrLevel::FillPatch(*this, Sborder, 4, prev_time, State_Type, 0, NUM_STATE);
    clean_state(Sborder);

    // Initialize the new-time data from the old time data.

    MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

    // Do the hydro update.

    // Construct the primitive variables.
    cons_to_prim(time);

    // Construct the hydro source.
    construct_hydro_source(time, dt);
    MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NUM_STATE, 0);
    clean_state(S_new);

    // Clean up our temporary MultiFab.
    Sborder.clear();

    return dt;

}
