
#include <iostream>
#include <Castro.H>

amrex::LevelBld* getLevelBld ();

int
main (int argc, char* argv[])
{
    // Disable AMReX verbosity for the purposes of this benchmark.
    amrex::SetVerbose(0);

    // Initialize AMReX.
    amrex::Initialize(argc,argv);

    // Collect n_cell now. We're only actually going to run
    // the code if n_cell > 0. If the user doesn't set n_cell,
    // we'll just print out some information about the code.
    amrex::ParmParse pp;

    int n_cell = 0;
    pp.query("n_cell", n_cell);
    
    if (n_cell <= 0)
    {
        amrex::Print() << std::endl;
        amrex::Print() << "mini-Castro is an astrophysical 3D Sedov-Taylor blast wave benchmark." << std::endl;
        amrex::Print() << "The algorithm closely resembles that of the Castro code." << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "This mini-app has several inputs parameters (values in parentheses are the defaults):" << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "n_cell: The number of zones per dimension. This parameter must be set to run the code." << std::endl;
        amrex::Print() << "max_box_size (64): The maximum size of a box in the domain." << std::endl;
        amrex::Print() << "min_box_size (16): The minimum size of a box in the domain." << std::endl;
        amrex::Print() << "max_level (0): The maximum adaptive mesh refinement level (zero-indexed)." << std::endl;
        amrex::Print() << "stop_time (0.01): The stopping time of the simulation, in seconds." << std::endl;
        amrex::Print() << "max_step (10000000): The maximum number of timesteps to take." << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "Example program invocation:" << std::endl;
        amrex::Print() << "./mini-Castro3d.pgi.MPI.CUDA.ex n_cell = 128 max_box_size = 128 min_box_size = 32 max_level = 0" << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "n_cell should be used to control the total amount of work to do. There are n_cell^3" << std::endl <<
                          "zones in the simulation, subdivided into boxes of various sizes. The minimum size" << std::endl <<
                          "of a box in any dimension is set by min_box_size, and similarly for max_box_size." << std::endl <<
                          "For appropriate load balancing, there must be enough boxes as there are MPI ranks" << std::endl <<
                          "so that every rank has work to do. Ideally there are multiple boxes per rank." << std::endl <<
                          "For example, if n_cell = 768, and you have 128 MPI ranks, then a reasonable choice" << std::endl <<
                          "would be max_box_size = 96, since (768/96)^3 = 512 = 128 * 4, so that every rank" << std::endl <<
                          "has 4 boxes to work with. min_box_size can be increased to decrease the number of boxes." << std::endl <<
                          "Generally speaking, larger boxes are more efficient to compute on. Appropriate choices" << std::endl <<
                          "for min_box_size and max_box_size usually must be empirically determined by trying" << std::endl <<
                          "multiple values and finding the fastest run time." << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "The simulation prints a Figure of Merit at the end which measures the simulation throughput." << std::endl;
        amrex::Print() << "The FOM measures the average number of zones advanced per microsecond (higher is better)." << std::endl;
        amrex::Print() << "To disable printing the FOM, set fom = 0." << std::endl;
        amrex::Print() << std::endl;
        amrex::Print() << "To track the state of the simulation, the effective radius of the blast wave is periodically calculated and printed." << std::endl;
        amrex::Print() << std::endl;
    }
    else
    {
        amrex::Print() << std::endl << "Starting mini-Castro..." << std::endl << std::endl;

        int max_step = 10000000;
        amrex::Real stop_time = 1.0e-2;
        int do_fom = 1;

        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
        pp.query("fom", do_fom);

        // Set the geometry parameters for this problem.
        // They are hardcoded for the Sedov blast wave
        // that we are solving.

        amrex::ParmParse pp_geom("geometry");

        std::vector<int> periodic{1, 1, 1};
        std::vector<amrex::Real> prob_lo{0.0, 0.0, 0.0};
        std::vector<amrex::Real> prob_hi{1.0e9, 1.0e9, 1.0e9};

        pp_geom.add("coord_sys", 0);
        pp_geom.addarr("is_periodic", periodic);
        pp_geom.addarr("prob_lo", prob_lo);
        pp_geom.addarr("prob_hi", prob_hi);

        // Use n_cell to replace amr.n_cell for a friendlier
        // user experience for those unfamiliar with AMReX.

        amrex::ParmParse pp_amr("amr");

        std::vector<int> n_cell_arr{n_cell, n_cell, n_cell};
        pp_amr.addarr("n_cell", n_cell_arr);

        // Use max_box_size to replace amr.max_grid_size.

        int max_box_size = 64;
        pp.query("max_box_size", max_box_size);
        pp_amr.add("max_grid_size", max_box_size);

        // Use min_box_size to replace amr.blocking_factor.

        int min_box_size = 16;
        pp.query("min_box_size", min_box_size);
        pp_amr.add("blocking_factor", min_box_size);

        // Use max_level to replace amr.max_level.

        int max_level = 0;
        pp.query("max_level", max_level);
        pp_amr.add("max_level", max_level);

        amrex::Print() << "Initializing AMR driver using the following runtime parameters:" << std::endl << std::endl;
        amrex::Print() << "n_cell = " << n_cell << std::endl;
        amrex::Print() << "max_box_size = " << max_box_size << std::endl;
        amrex::Print() << "min_box_size = " << min_box_size << std::endl;
        amrex::Print() << "max_level = " << max_level << std::endl;
        amrex::Print() << "max_step = " << max_step << std::endl;
        amrex::Print() << "stop_time = " << stop_time << std::endl;
        amrex::Print() << std::endl;

        amrex::Amr* amrptr = new amrex::Amr(getLevelBld());

        amrex::Print() << "Running simulation and calculating diagnostics every " << Castro::diagnostic_interval << " timesteps..." << std::endl << std::endl;

        amrptr->init(0.0, stop_time);

        amrex::Real dRunTime1 = amrex::ParallelDescriptor::second();

        while ( amrptr->okToContinue()                            &&
               (amrptr->levelSteps(0) < max_step || max_step < 0) &&
               (amrptr->cumTime() < stop_time || stop_time < 0.0) )
        {
            //
            // Do a timestep.
            //
            amrptr->coarseTimeStep(stop_time);

        }

        int nsteps = amrptr->levelSteps(0);

        // Start calculating the figure of merit for this run: average number of zones
        // advanced per microsecond. This must be done before we delete the Amr
        // object because we need to scale it by the number of zones on the coarse grid.

        long numPtsCoarseGrid = amrptr->getLevel(0).boxArray().numPts();
        amrex::Real fom = Castro::num_zones_advanced * numPtsCoarseGrid;

        delete amrptr;

        amrex::Real dRunTime2 = amrex::ParallelDescriptor::second();

        amrex::Real runtime = dRunTime2 - dRunTime1;

        const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        fom = fom / runtime / 1.e6;

        amrex::Print() << std::endl;
        amrex::Print() << "Simulation completed!" << std::endl;
        amrex::Print() << "Number of timesteps taken: " << nsteps << std::endl;
        amrex::Print() << std::endl;
        if (do_fom) {
            amrex::Print() << "Figure of Merit (zones / usec): " << std::fixed << std::setprecision(3) << fom << "\n";
            amrex::Print() << std::endl;
        }

    }

    amrex::Finalize();

    return 0;
}
