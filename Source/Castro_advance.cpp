
#include <Castro.H>
#include <Castro_F.H>

#include <cmath>
#include <climits>

using std::string;
using namespace amrex;

Real
Castro::advance (Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::advance()");

    Real dt_new = dt;

    // Prepare for the advance.

    initialize_advance(time, dt);

    // Do the advance from t = time to t = time + dt.

    for (int iter = 0; iter < MOL_STAGES; ++iter) {

        // Perform initialization steps.

        initialize_do_advance(time, dt, iter, MOL_STAGES);

        // Do the actual advance.

	dt_new = do_advance(time, dt, iter, MOL_STAGES);

        // Clean up.

        finalize_do_advance(time, dt, iter, MOL_STAGES);

    }

    // Clean up.

    finalize_advance(time, dt);

    return dt_new;
}

Real
Castro::do_advance (Real time, Real dt, int sub_iteration, int sub_ncycle)
{

    // This routine will advance the old state data (called S_old here)
    // to the new time, for a single level.  The new data is called
    // S_new here.

    BL_PROFILE("Castro::do_advance()");

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Check for NaN's.

//    check_for_nan(S_old);

    if (sub_iteration == 0) {

      // Initialize the new-time data.

      MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

      // Store the result of the burn in Sburn for later stages.
      // In the full Castro application, there would have been a
      // nuclear reactions burn for dt / 2 prior to this step.

      MultiFab::Copy(Sburn, Sborder, 0, 0, NUM_STATE, Sburn.nGrow());

      // We'll add each stage's contribution to -div{F(U)} as we compute them
      hydro_source.setVal(0.0);

    }

    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn, if it had occurred.

    construct_mol_hydro_source(time, dt, sub_iteration, sub_ncycle);

    // For MOL integration, we are done with this stage, unless it is
    // the last stage.

    if (sub_iteration == sub_ncycle-1) {

	// We just finished the last stage of the MOL integration.
	// Construct S_new now using the weighted sum of the updates,
	// starting from the post-burn state.

	MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
	MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, S_new.nComp(), 0);
      
	// Define the temperature now.
	clean_state(S_new);

	// Check for NaN's
//      check_for_nan(S_new);

    }

    return dt;

}



void
Castro::initialize_do_advance(Real time, Real dt, int sub_iteration, int sub_ncycle)
{
    BL_PROFILE("Castro::initialize_do_advance()");

    int finest_level = parent->finestLevel();

    // For the hydrodynamics update we need to have 4 ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    Sborder.define(grids, dmap, NUM_STATE, 4);

    Real expand_time;

    if (sub_iteration == 0) {

	// First MOL stage
	expand_time = state[State_Type].prevTime();

    } else {

	// The initial state for the kth stage follows the Butcher
	// tableau.  We need to create the proper state starting with
	// the result after the first nuclear reactions burn (which we
        // copied into Sburn) and we need to fill ghost cells. As above,
        // we are not actually doing a burn in this mini-app, but we
        // do this step for completeness. We'll overwrite S_new with this
        // information, since we don't need it anymore.
	MultiFab& S_new = get_new_data(State_Type);

        std::vector< std::vector<Real> > a_mol{ {0.0, 0.0}, {1.0, 0.0} };

	MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
	for (int i = 0; i < sub_iteration; ++i)
	    MultiFab::Saxpy(S_new, dt*a_mol[sub_iteration][i], *k_mol[i], 0, 0, S_new.nComp(), 0);

	expand_time = state[State_Type].curTime();

    }

    expand_state(Sborder, expand_time, 4);
}



void
Castro::finalize_do_advance(Real time, Real dt, int sub_iteration, int sub_ncycle)
{
    BL_PROFILE("Castro::finalize_do_advance()");

    Sborder.clear();
}



void
Castro::initialize_advance(Real time, Real dt)
{
    BL_PROFILE("Castro::initialize_advance()");

    // Swap the new data from the last timestep into the old state data.

    state[State_Type].allocOldData();
    state[State_Type].swapTimeLevels(dt);

    // Ensure data is valid before beginning advance.

    clean_state(get_old_data(State_Type));

    // This array holds the hydrodynamics update.

    hydro_source.define(grids,dmap,NUM_STATE,0);

    k_mol.resize(MOL_STAGES);
    for (int n = 0; n < MOL_STAGES; ++n) {
	k_mol[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
	k_mol[n]->setVal(0.0);
    }

    // This array holds the state that would be updated after
    // performing a nuclear reactions step. In mini-Castro we
    // do not do nuclear reactions for simplicity, but we
    // maintain the code flow for representativeness.

    Sburn.define(grids, dmap, NUM_STATE, 0);

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	fluxes[dir]->setVal(0.0);
}



void
Castro::finalize_advance(Real time, Real dt)
{
    BL_PROFILE("Castro::finalize_advance()");

    // Update flux registers.

    if (level < parent->finestLevel()) {
        for (int i = 0; i < 3; ++i) {
            getLevel(level+1).flux_reg.CrseInit(*fluxes[i], i, 0, 0, NUM_STATE, -1.0);
        }
    }

    if (level > 0) {
        for (int i = 0; i < 3; ++i) {
            flux_reg.FineAdd(*fluxes[i], i, 0, 0, NUM_STATE, 1.0);
        }
    }

    hydro_source.clear();
    k_mol.clear();
    Sburn.clear();

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts();
}
