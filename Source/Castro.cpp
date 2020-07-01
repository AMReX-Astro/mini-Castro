
#include <vector>
#include <iostream>
#include <string>

#include <AMReX_Utility.H>
#include <Castro.H>
#include <Castro_F.H>
#include <AMReX_TagBox.H>
#include <AMReX_FillPatchUtil.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

using namespace amrex;

Real Castro::num_zones_advanced = 0.0;
int Castro::diagnostic_interval = 50;

// Choose tile size based on whether we're using a GPU.

#ifdef AMREX_USE_GPU
IntVect Castro::tile_size(1024000, 1024000, 1024000);
#else
IntVect Castro::tile_size(1024, 16, 16);
#endif

void
Castro::variableCleanUp ()
{
    desc_lst.clear();

    eos_finalize();
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

    if (lev > 0) {
        amrex::Abort("AMR is unsupported");
    }

    // Initialize volume and area arrays.

    volume.clear();
    volume.define(grids, dmap, 1, 4);

    auto dx = geom.CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(volume, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        auto vol = volume[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            set_volume(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                       AMREX_ARR4_TO_FORTRAN_ANYD(vol),
                       AMREX_ZFILL(dx.data()));
        });

    }

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
	area[dir].define(getEdgeBoxArray(dir), dmap, 1, 4);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(area[dir], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            auto ar = area[dir][mfi].array();

            CASTRO_LAUNCH_LAMBDA(box, lbx,
            {
                set_area(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                         AMREX_ARR4_TO_FORTRAN_ANYD(ar),
                         AMREX_ZFILL(dx.data()),
                         dir + 1);
            });
        }
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

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        auto state_arr = S_new[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            initdata(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
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

    const auto dx = geom.CellSizeArray();

    Real* dt_loc = static_cast<Real*>(amrex::The_Managed_Arena()->alloc(sizeof(Real)));

    *dt_loc = std::numeric_limits<amrex::Real>::max();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(stateMF, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();

        auto state_arr = stateMF[mfi].array();

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            estdt(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                  AMREX_ARR4_TO_FORTRAN_ANYD(state_arr),
                  AMREX_ZFILL(dx.data()), dt_loc);
        });
    }

    Real dt = *dt_loc;

    amrex::The_Managed_Arena()->free(dt_loc);

    // Reduce over all MPI ranks.
    ParallelDescriptor::ReduceRealMin(dt);

    const Real cfl = 0.5;
    dt *= cfl;

    return dt;
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

    Real change_max = 1.1e0;

    Real dt_0 = 1.0e+100;
    dt_min[0] = estTimeStep(dt_level[0]);

    //
    // Limit dt by change_max * old dt
    //
    dt_min[0] = std::min(dt_min[0], change_max * dt_level[0]);
    dt_0 = std::min(dt_0, dt_min[0]);

    //
    // Limit dt by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    dt_level[0] = dt_0;
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

    Real dt_0 = 1.0e+100;
    dt_level[0] = getLevel(0).initialTimeStep();
    dt_0 = std::min(dt_0, dt_level[0]);

    //
    // Limit dt by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    dt_level[0] = dt_0;
}

void
Castro::post_timestep (int iteration)
{
    BL_PROFILE("Castro::post_timestep()");

    if (parent->levelSteps(0) % diagnostic_interval == 0)
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

        // Get device pointers to the reduction variables.
        Real* blast_mass_loc = static_cast<Real*>(amrex::The_Managed_Arena()->alloc(sizeof(Real)));
        Real* blast_radius_loc = static_cast<Real*>(amrex::The_Managed_Arena()->alloc(sizeof(Real)));

        *blast_mass_loc = 0.0;
        *blast_radius_loc = 0.0;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, tile_size); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();

            auto state_arr = S_new[mfi].array();

            CASTRO_LAUNCH_LAMBDA(box, lbx,
            {
                calculate_blast_radius(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                       AMREX_ARR4_TO_FORTRAN_ANYD(state_arr),
                                       AMREX_ZFILL(dx.data()), AMREX_ZFILL(problo.data()), AMREX_ZFILL(probhi.data()),
                                       blast_mass_loc, blast_radius_loc, max_density);
            });

        }

        Real blast_mass = *blast_mass_loc;
        Real blast_radius = *blast_radius_loc;

        amrex::The_Managed_Arena()->free(blast_mass_loc);
        amrex::The_Managed_Arena()->free(blast_radius_loc);

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
}

void
Castro::post_init (Real stop_time)
{}

void
Castro::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{}

// Given State_Type state data, perform a number of cleaning steps to make
// sure the data is sensible.

void
Castro::clean_state(MultiFab& state)
{
    BL_PROFILE("Castro::clean_state()");

    int ng = state.nGrow();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, tile_size); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.growntilebox(ng);
        auto state_arr = state[mfi].array();

        // Ensure the density is larger than the density floor.

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            enforce_minimum_density(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                                    AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Ensure all species are normalized.

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            normalize_species(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                              AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Ensure (rho e) isn't too small or negative

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            reset_internal_e(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                             AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });

        // Make the temperature be consistent with the internal energy.

        CASTRO_LAUNCH_LAMBDA(box, lbx,
        {
            compute_temp(AMREX_ARLIM_ANYD(lbx.loVect()), AMREX_ARLIM_ANYD(lbx.hiVect()),
                         AMREX_ARR4_TO_FORTRAN_ANYD(state_arr));
        });
    }
}
