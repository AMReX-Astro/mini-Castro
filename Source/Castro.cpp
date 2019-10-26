
#include <vector>
#include <iostream>
#include <string>

#include <AMReX_Utility.H>
#include <Castro.H>
#include <Castro_F.H>
#include <AMReX_TagBox.H>
#include <AMReX_FillPatchUtil.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

long Castro::num_zones_advanced = 0;
int Castro::diagnostic_interval = 50;

void
Castro::variableCleanUp ()
{
    desc_lst.clear();

    eos_finalize();

    network_finalize();

    probinit_finalize();
}

Castro::Castro (Amr&            papa,
                int             lev,
                const Geometry& level_geom,
                const BoxArray& bl,
		const DistributionMapping& dm,
                Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{

    BL_PROFILE("Castro::Castro()");

    // Initialize volume, area, flux arrays.

    volume.clear();
    volume.define(grids, dmap, 1, 4);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
	area[dir].define(getEdgeBoxArray(dir), dmap, 1, 4);
        geom.GetFaceArea(area[dir],dir);
    }

    fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, NUM_STATE, 0));

    if (level > 0) {

	flux_reg.define(grids, dmap, crse_ratio, level, NUM_STATE);
	flux_reg.setVal(0.0);

    }

}

Castro::~Castro ()
{
}

void
Castro::initData ()
{
    BL_PROFILE("Castro::initData()");

    // Loop over grids and initialize data.

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    S_new.setVal(0.);

    // make sure dx = dy = dz -- that's all we guarantee to support
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
    {
	amrex::Abort("We don't support dx != dy != dz");
    }

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        auto state_arr = S_new[mfi].array();

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_initdata
                (AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                 AMREX_ARR4_TO_FORTRAN_ANYD(state_arr), AMREX_ZFILL(dx.data()),
                 AMREX_ZFILL(problo.data()), AMREX_ZFILL(probhi.data()));
        });
    }
}

void
Castro::init (AmrLevel &old)
{
    BL_PROFILE("Castro::init(old)");

    Castro* oldlev = (Castro*) &old;

    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& state_MF = get_new_data(State_Type);
    FillPatch(old, state_MF, state_MF.nGrow(), cur_time, State_Type, 0, state_MF.nComp());
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Castro::init ()
{
    BL_PROFILE("Castro::init()");

    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    Real time = cur_time;

    setTimeLevel(time, dt_old, dt);

    MultiFab& state_MF = get_new_data(State_Type);
    FillCoarsePatch(state_MF, 0, cur_time, State_Type, 0, state_MF.nComp());
}

Real
Castro::initialTimeStep ()
{
    BL_PROFILE("Castro::initialTimeStep()");

    const Real dummy_dt = 0.0;
    const Real init_shrink = 0.01;

    Real init_dt = init_shrink * estTimeStep(dummy_dt);

    return init_dt;
}

Real
Castro::estTimeStep (Real dt_old)
{
    BL_PROFILE("Castro::estTimeStep()");

    const MultiFab& stateMF = get_new_data(State_Type);

    amrex::Real estdt = std::numeric_limits<amrex::Real>::max();

    const auto dx = geom.CellSizeArray();

    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();

        auto state_arr = stateMF[mfi].array();

        // Get a device pointer for estdt.
        Real* estdt_loc = AMREX_MFITER_REDUCE_MIN(&estdt);

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_estdt(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                     AMREX_ARR4_TO_FORTRAN_ANYD(state_arr),
                     AMREX_ZFILL(dx.data()), estdt_loc);
        });
    }

    // Reduce over all MPI ranks.
    ParallelDescriptor::ReduceRealMin(estdt);

    const Real cfl = 0.5;
    estdt *= cfl;

    return estdt;
}

void
Castro::computeNewDt (int                   finest_level,
                      int                   sub_cycle,
                      Vector<int>&           n_cycle,
                      const Vector<IntVect>& ref_ratio,
                      Vector<Real>&          dt_min,
                      Vector<Real>&          dt_level,
                      Real                  stop_time,
                      int                   post_regrid_flag)
{
    BL_PROFILE("Castro::computeNewDt()");

    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    Real change_max = 1.1e0;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Castro& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1)
    {
       //
       // Limit dt's by pre-regrid dt
       //
       for (i = 0; i <= finest_level; i++)
       {
           dt_min[i] = std::min(dt_min[i],dt_level[i]);
       }
    }
    else
    {
       //
       // Limit dt's by change_max * old dt
       //
       for (i = 0; i <= finest_level; i++)
       {
           dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
       }
    }

    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::computeInitialDt (int                   finest_level,
                          int                   sub_cycle,
                          Vector<int>&           n_cycle,
                          const Vector<IntVect>& ref_ratio,
                          Vector<Real>&          dt_level,
                          Real                  stop_time)
{
    BL_PROFILE("Castro::computeInitialDt()");

    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    ///TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::post_timestep (int iteration)
{
    BL_PROFILE("Castro::post_timestep()");

    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    // Now do the refluxing.

    if (level < parent->finestLevel())
	reflux(level, level+1);

    // Ensure consistency with finer grids.

    if (level < finest_level)
	avgDown();

    MultiFab& S_new = get_new_data(State_Type);

    // Clean up any aberrant state data generated by the reflux and average-down,
    // and then update quantities like temperature to be consistent.

    clean_state(S_new);

    if (level == 0 && parent->levelSteps(0) % diagnostic_interval == 0)
    {
        // As a diagnostic quantity, we'll print the current blast radius
        // (the location of the shock). We can estimate this using the
        // location where the density on the domain is maximum. Then we
        // can refine the estimate by doing a weighted sum on the domain.

        // We'll assume that everything we care about is on the finest level.
        MultiFab& S_new = getLevel(parent->finestLevel()).get_new_data(State_Type);
        Real max_density = S_new.max(Density);

        const auto dx = getLevel(parent->finestLevel()).geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();
        const auto probhi = geom.ProbHiArray();

        Real blast_radius = 0.0;
        Real blast_mass = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(sum:blast_mass,blast_radius)
#endif
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();

            auto state_arr = S_new[mfi].array();

            // Get device pointers to the reduction variables.
            Real* blast_mass_loc = AMREX_MFITER_REDUCE_SUM(&blast_mass);
            Real* blast_radius_loc = AMREX_MFITER_REDUCE_SUM(&blast_radius);

            AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
            {
                calculate_blast_radius(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                       AMREX_ARR4_TO_FORTRAN_ANYD(state_arr),
                                       AMREX_ZFILL(dx.data()), AMREX_ZFILL(problo.data()), AMREX_ZFILL(probhi.data()),
                                       blast_mass_loc, blast_radius_loc, max_density);
            });

        }

        // Reduce over MPI ranks.
        amrex::ParallelDescriptor::ReduceRealSum(blast_mass);
        amrex::ParallelDescriptor::ReduceRealSum(blast_radius);

        amrex::Print() << std::scientific << std::setprecision(6) << "Blast radius at step " << parent->levelSteps(0) << ", time " << state[State_Type].curTime()
                       << ": " << std::fixed << std::setprecision(3) << (blast_radius / blast_mass) / 1.0e5 << " km" << std::endl;
    }

}

void
Castro::post_regrid (int lbase, int new_finest)
{
    BL_PROFILE("Castro::post_regrid()");
}

void
Castro::post_init (Real stop_time)
{
    BL_PROFILE("Castro::post_init()");

    if (level > 0) return;

    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}



void
Castro::reflux(int crse_level, int fine_level)
{
    BL_PROFILE("Castro::reflux()");

    BL_ASSERT(fine_level > crse_level);

    const Real strt = ParallelDescriptor::second();

    FluxRegister* reg;

    for (int lev = fine_level; lev > crse_level; --lev) {

	reg = &getLevel(lev).flux_reg;

	Castro& crse_lev = getLevel(lev-1);
	Castro& fine_lev = getLevel(lev);

	MultiFab& state = crse_lev.get_new_data(State_Type);

	// Clear out the data that's not on coarse-fine boundaries so that this register only
	// modifies the fluxes on coarse-fine interfaces.

	reg->ClearInternalBorders(crse_lev.geom);

	// Trigger the actual reflux on the coarse level now.

	reg->Reflux(state, crse_lev.volume, 1.0, 0, 0, NUM_STATE, crse_lev.geom);

	// We no longer need the flux register data, so clear it out.

	reg->setVal(0.0);

    }

}

void
Castro::avgDown ()
{
    BL_PROFILE("Castro::avgDown()");

    if (level == parent->finestLevel()) return;

    Castro& fine_lev = getLevel(level+1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom =          geom;

    MultiFab& S_crse = get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse,
                        fgeom, cgeom,
                        0, S_fine.nComp(), fine_ratio);
}

void
Castro::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{
    BL_PROFILE("Castro::errorEst()");

    auto mf = derive("density", time, 1);

    BL_ASSERT(mf);

    const int8_t set   = (int8_t) TagBox::SET;
    const int8_t clear = (int8_t) TagBox::CLEAR;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        auto tags_arr = tags[mfi].array();
        auto data_arr = (*mf)[mfi].array();

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_denerror(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                        (int8_t*) AMREX_ARR4_TO_FORTRAN_ANYD(tags_arr),
                        AMREX_ARR4_TO_FORTRAN_ANYD(data_arr),
                        set, clear);
        });
    }
}

// Fill a version of the state with ng ghost zones from the state data.
void
Castro::expand_state(MultiFab& S, Real time, int ng)
{
    BL_PROFILE("Castro::expand_state()");

    BL_ASSERT(S.nGrow() >= ng);

    AmrLevel::FillPatch(*this,S,ng,time,State_Type,0,NUM_STATE);

    clean_state(S);
}


void
Castro::check_for_nan(MultiFab& state, int check_ghost)
{
  BL_PROFILE("Castro::check_for_nan()");

  int ng = 0;
  if (check_ghost == 1) {
    ng = state.nComp();
  }

  if (state.contains_nan(Density,state.nComp(),ng,true))
    {
      for (int i = 0; i < state.nComp(); i++)
        {
	  if (state.contains_nan(Density + i, 1, ng, true))
            {
	      std::string abort_string = std::string("State has NaNs in the ") + desc_lst[State_Type].name(i) + std::string(" component::check_for_nan()");
	      amrex::Abort(abort_string.c_str());
            }
        }
    }
}

// Given State_Type state data, perform a number of cleaning steps to make
// sure the data is sensible.

void
Castro::clean_state(MultiFab& state)
{
    BL_PROFILE("Castro::clean_state()");

    int ng = state.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.growntilebox(ng);
        auto state_arr = state[mfi].array();

        // Ensure the density is larger than the density floor.

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_enforce_minimum_density(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                       AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Ensure all species are normalized.

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_normalize_species(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                 AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Ensure (rho e) isn't too small or negative

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_reset_internal_e(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Make the temperature be consistent with the internal energy.

        AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx,
        {
            ca_compute_temp(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                            AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });
    }
}
