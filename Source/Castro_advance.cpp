
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
	dt_new = do_advance(time, dt, iter, MOL_STAGES);
    }

    // Clean up.

    finalize_advance(time, dt);

    return dt_new;
}

Real
Castro::do_advance (Real time, Real dt, int sub_iteration, int sub_ncycle)
{

  // this routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.

    BL_PROFILE("Castro::do_advance()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Perform initialization steps.

    initialize_do_advance(time, dt, sub_iteration, sub_ncycle);

    // Check for NaN's.

//    check_for_nan(S_old);

    if (sub_iteration == 0) {

      // Initialize the new-time data.

      MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

      // store the result of the burn in Sburn for later stages

      MultiFab::Copy(Sburn, Sborder, 0, 0, NUM_STATE, Sburn.nGrow());

      // we'll add each stage's contribution to -div{F(U)} as we compute them
      hydro_source.setVal(0.0);

    }

    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn 

    construct_mol_hydro_source(time, dt, sub_iteration, sub_ncycle);

    // For MOL integration, we are done with this stage, unless it is
    // the last stage

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

    finalize_do_advance(time, dt, sub_iteration, sub_ncycle);

    return dt;

}



void
Castro::initialize_do_advance(Real time, Real dt, int sub_iteration, int sub_ncycle)
{

    BL_PROFILE_VAR("Castro::initialize_do_advance()", CA_INIT_DO_ADV);

    int finest_level = parent->finestLevel();

    // For the hydrodynamics update we need to have 4 ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    if (sub_iteration == 0) {

	// first MOL stage
	Sborder.define(grids, dmap, NUM_STATE, 4);
	const Real prev_time = state[State_Type].prevTime();
	expand_state(Sborder, prev_time, 4);

    } else {

	// the initial state for the kth stage follows the Butcher
	// tableau.  We need to create the proper state starting with
	// the result after the first dt/2 burn (which we copied into
	// Sburn) and we need to fill ghost cells.  

	// We'll overwrite S_old with this information, since we don't
	// need it anymorebuild this state temporarily in S_new (which
	// is State_Data) to allow for ghost filling.
	MultiFab& S_new = get_new_data(State_Type);

	MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
	for (int i = 0; i < sub_iteration; ++i)
	    MultiFab::Saxpy(S_new, dt*a_mol[sub_iteration][i], *k_mol[i], 0, 0, S_new.nComp(), 0);

	Sborder.define(grids, dmap, NUM_STATE, 4);
	const Real new_time = state[State_Type].curTime();
	expand_state(Sborder, new_time, 4);

    }

    BL_PROFILE_VAR_STOP(CA_INIT_DO_ADV);

}



void
Castro::finalize_do_advance(Real time, Real dt, int sub_iteration, int sub_ncycle)
{
    BL_PROFILE_VAR("Castro::finalize_do_advance()", CA_FIN_DO_ADV);

    Sborder.clear();

    BL_PROFILE_VAR_STOP(CA_FIN_DO_ADV);
}



void
Castro::initialize_advance(Real time, Real dt)
{
    BL_PROFILE_VAR("Castro::initialize_advance()", CA_INIT_ADV);

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

    // for the post-burn state
    Sburn.define(grids, dmap, NUM_STATE, 0);

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	fluxes[dir]->setVal(0.0);

    BL_PROFILE_VAR_STOP(CA_INIT_ADV);
}



void
Castro::finalize_advance(Real time, Real dt)
{
    BL_PROFILE_VAR("Castro::finalize_advance()", CA_FIN_ADV);

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

    BL_PROFILE_VAR_STOP(CA_FIN_ADV);
}
