
#ifndef _Castro_H_
#define _Castro_H_

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FluxRegister.H>

#define NumSpec 13
#define MOL_STAGES 2
#define NQAUX 4
#define NGDNV 6
#define QVAR 21

enum StateType { State_Type };

enum Conserved { Density = 0, Xmom, Ymom, Zmom, Eden, Eint, Temp, FirstSpec, NUM_STATE = FirstSpec + NumSpec };

#define AMREX_ARR4_TO_FORTRAN_ANYD(a) a.p,&((a).begin.x),amrex::GpuArray<int,3>{(a).end.x-1,(a).end.y-1,(a).end.z-1}.data()

// If we are using OpenACC, disable AMReX from launching CUDA kernels for our code.

#ifdef AMREX_USE_ACC
#define CASTRO_LAUNCH_LAMBDA(box, lbx, lambda) { Box lbx = box; lambda; }
#define CASTRO_MFITER_REDUCE_MIN(x) x
#define CASTRO_MFITER_REDUCE_SUM(x) x
#else
#define CASTRO_LAUNCH_LAMBDA(box, lbx, lambda) AMREX_LAUNCH_DEVICE_LAMBDA(box, lbx, lambda)
#define CASTRO_MFITER_REDUCE_MIN(x) AMREX_MFITER_REDUCE_MIN(x)
#define CASTRO_MFITER_REDUCE_SUM(x) AMREX_MFITER_REDUCE_SUM(x)
#endif

class Castro
    :
    public amrex::AmrLevel
{
public:

    // Default constructor; builds invalid object
    Castro () {}

    // The basic constructor
    Castro (amrex::Amr&                       papa,
            int                               lev,
            const amrex::Geometry&            level_geom,
            const amrex::BoxArray&            bl,
            const amrex::DistributionMapping& dm,
            amrex::Real                       time);

    // The destructor
    virtual ~Castro () override;

    // Get a reference to a given level
    Castro& getLevel (int lev);

    // Define data descriptors
    static void variableSetUp ();

    // Cleanup data descriptors at end of run
    static void variableCleanUp ();

    // Initialize grid data at problem start-up
    virtual void initData () override;

    // Initialize data on this level from another Castro (during regrid)
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;

    // Advance grids at this level in time
    virtual amrex::Real advance (amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

    // Convert conservative variables to primitive variables
    void cons_to_prim(const amrex::Real time);
    
    // Construct the hydrodynamic source term
    void construct_hydro_source(amrex::Real time, amrex::Real dt);

    // Estimate time step
    amrex::Real estTimeStep (amrex::Real dt_old);

    // Compute initial time step
    amrex::Real initialTimeStep ();

    // Compute initial timestep
    virtual void computeInitialDt (int                                  finest_level,
                                   int                                  sub_cycle,
                                   amrex::Vector<int>&                  n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                          stop_time) override;

    // Compute new timestep
    virtual void computeNewDt (int                                  finest_level,
                               int                                  sub_cycle,
                               amrex::Vector<int>&                  n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                          stop_time,
                               int                                  post_regrid_flag) override;

    // Do work after timestep()
    virtual void post_timestep (int iteration) override;

    // Do work after regrid()
    virtual void post_regrid (int lbase, int new_finest) override;   

    // Do work after init()
    virtual void post_init (amrex::Real stop_time) override;

    // Error estimation for regridding
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
			   int                 n_error_buf = 0,
			   int                 ngrow = 0) override;

    // Apply a number of corrections to ensure consistency in the state
    void clean_state (amrex::MultiFab& state);

    // Initialize the flux register with coarse data.
    void FluxRegCrseInit();

    // Update the flux register with fine data.
    void FluxRegFineAdd();
    
    // Update coarse levels with flux correction from fine levels
    void reflux (int crse_level, int fine_level);

    // Average down data from fine levels onto underlying coarse levels
    void avgDown ();

    // A record of how many cells we have advanced throughout the simulation.
    // This is saved as a real because we will be storing the number of zones
    // advanced as a ratio with the number of zones on the coarse grid (to
    // help prevent the number from getting too large), and that may be a
    // non-integer number.
    static amrex::Real num_zones_advanced;

    // How often should we print out diagnostic output?
    static int diagnostic_interval;

protected:

    // A state array with ghost zones
    amrex::MultiFab Sborder;

    // The primitive variable state array.
    amrex::MultiFab q;

    // The auxiliary primitive variable state array.
    amrex::MultiFab qaux;

    // Source term representing hydrodynamics update
    amrex::MultiFab hydro_source;

    // Hydrodynamic fluxes
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
    amrex::FluxRegister flux_reg;

    // Geometric data
    amrex::MultiFab volume;
    amrex::MultiFab area[3];

    static amrex::IntVect tile_size;

};

inline
Castro&
Castro::getLevel (int lev)
{
    return *(Castro *) &parent->getLevel(lev);
}

#endif /*_Castro_H_*/
