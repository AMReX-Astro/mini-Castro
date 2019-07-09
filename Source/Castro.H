
#ifndef _Castro_H_
#define _Castro_H_

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#include <memory>
#include <iostream>

using std::istream;
using std::ostream;

#define NumSpec 13

enum StateType { State_Type };

enum Conserved { Density = 0, Xmom, Ymom, Zmom, Eden, Eint, Temp, FirstSpec, NUM_STATE = FirstSpec + NumSpec };

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar media
//

class Castro
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Castro () {}
    //
    //The basic constructor.
    //
    Castro (amrex::Amr&            papa,
            int                    lev,
            const amrex::Geometry& level_geom,
            const amrex::BoxArray& bl,
            const amrex::DistributionMapping& dm,
            amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~Castro () override;
    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Define tagging functions.
    //
    static void ErrorSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Initialize grid data at problem start-up.
    //
    virtual void initData () override;
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (amrex::Real time,
                               amrex::Real dt_old,
                               amrex::Real dt_new) override;
    //
    // Initialize data on this level from another Castro (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    // Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    // Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                          amrex::Real dt,
                          int  iteration,
                          int  ncycle) override;

    amrex::Real do_advance (amrex::Real time,
                     amrex::Real dt,
                     int  amr_iteration,
                     int  amr_ncycle,
                     int  sub_iteration,
                     int  sub_ncycle);

    amrex::Real retry_advance (amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

    void initialize_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

    void finalize_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

    void initialize_do_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void finalize_do_advance(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void construct_hydro_source(amrex::Real time, amrex::Real dt);

    void construct_mol_hydro_source(amrex::Real time, amrex::Real dt, int istage, int nstages);

    void check_for_nan(amrex::MultiFab& state, int check_ghost=0);

    //
    // Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    // Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    // Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    // Allocate data at old time.
    //
    virtual void allocOldData () override;
    //
    // Remove data at old time.
    //
    virtual void removeOldData () override;
    //
    // Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;
    //
    // Do work after regrid().
    //
    virtual void post_regrid (int lbase,
                              int new_finest) override;   
    //
    // Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;

    //
    // Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real  time,
			   int          n_error_buf = 0,
			   int          ngrow = 0) override;

    void reset_internal_energy (amrex::MultiFab& State);

    void computeTemp (amrex::MultiFab& State);

    void add_force_to_sources (amrex::MultiFab& force, amrex::MultiFab& sources, amrex::MultiFab& state);

    void apply_source_to_state (amrex::MultiFab& state, amrex::MultiFab& source, amrex::Real dt);

    void expand_state(amrex::MultiFab& S, amrex::Real time, int ng);

    static int QVAR;
    static int NQAUX;
    static int NQ;
    static int NGDNV;

    //
    // A record of how many cells we have advanced throughout the simulation
    //
    static long num_zones_advanced;

protected:

    //
    // Build a mask that ghost cells overlapping with interior cells in the same multifab
    // are set to 0, whereas others are set to 1.
    //
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > ib_mask;
    amrex::iMultiFab& build_interior_boundary_mask (int ng);

    //
    // A state array with ghost zones.
    //
    amrex::MultiFab Sborder;

    //
    // A state array for the post burn state to be used in MOL integration
    //
    amrex::MultiFab Sburn;

    //
    // Source term representing hydrodynamics update.
    //
    amrex::MultiFab hydro_source;

    //
    // Hydrodynamic (and radiation) fluxes.
    //
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
    amrex::FluxRegister flux_reg;

    // Scalings for the flux registers.
    amrex::Real flux_crse_scale;
    amrex::Real flux_fine_scale;

    //
    // Storage for the method of lines stages
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > k_mol;
    //

    // MOL Butcher tableau
    static int MOL_STAGES;
  
    static amrex::Vector< amrex::Vector<amrex::Real> > a_mol;
    static amrex::Vector<amrex::Real> b_mol;
    static amrex::Vector<amrex::Real> c_mol;

    static void network_init ();

    static void read_params ();

    Castro& getLevel (int lev);

    void FluxRegCrseInit();
    void FluxRegFineAdd();

    void reflux (int crse_level, int fine_level);

    void normalize_species (amrex::MultiFab& S_new);

    void enforce_min_density (amrex::MultiFab& S_old, amrex::MultiFab& S_new);

    void clean_state (amrex::MultiFab& state);

    void avgDown ();

    //
    // The data.
    //
    amrex::MultiFab volume;
    amrex::MultiFab area[3];

    //
    // Keep track of which AMR iteration we're on.
    //
    int iteration;

    //
    // Static data members.
    //
    static int NUM_GROW;

    static amrex::IntVect hydro_tile_size;

    static int num_state_type;

};

inline
Castro&
Castro::getLevel (int lev)
{
    return *(Castro *) &parent->getLevel(lev);
}

#endif /*_Castro_H_*/
