#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::sum_integrated_quantities ()
{

    if (verbose <= 0) return;

    bool local_flag = true;

    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real mom[3]      = { 0.0 };
    Real ang_mom[3]  = { 0.0 };
    Real com[3]      = { 0.0 };
    Real com_vel[3]  = { 0.0 };
    Real rho_e       = 0.0;
    Real rho_K       = 0.0;
    Real rho_E       = 0.0;
    int datwidth     = 14;
    int datprecision = 6;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass   += ca_lev.volWgtSum("density", time, local_flag);
        mom[0] += ca_lev.volWgtSum("xmom", time, local_flag);
	mom[1] += ca_lev.volWgtSum("ymom", time, local_flag);
	mom[2] += ca_lev.volWgtSum("zmom", time, local_flag);

	if (show_center_of_mass) {
	   com[0] += ca_lev.locWgtSum("density", time, 0, local_flag);
	   com[1] += ca_lev.locWgtSum("density", time, 1, local_flag);
	   com[2] += ca_lev.locWgtSum("density", time, 2, local_flag);
	}

       rho_e += ca_lev.volWgtSum("rho_e", time, local_flag);
       rho_E += ca_lev.volWgtSum("rho_E", time, local_flag);
    }
 
    if (verbose > 0)
    {

       const int nfoo = 10;
	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], ang_mom[0], ang_mom[1], ang_mom[2],
			  rho_e, rho_K, rho_E};
#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

	ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

	if (show_center_of_mass)
	    ParallelDescriptor::ReduceRealSum(com, 3, ParallelDescriptor::IOProcessorNumber());

	if (ParallelDescriptor::IOProcessor()) {

	    int i = 0;
	    mass       = foo[i++];
	    mom[0]     = foo[i++];
            mom[1]     = foo[i++];
            mom[2]     = foo[i++];
	    ang_mom[0] = foo[i++];
	    ang_mom[1] = foo[i++];
	    ang_mom[2] = foo[i++];
	    rho_e      = foo[i++];
	    rho_K      = foo[i++];
            rho_E      = foo[i++];
	    std::cout << '\n';
	    std::cout << "TIME= " << time << " MASS        = "   << mass      << '\n';
	    std::cout << "TIME= " << time << " XMOM        = "   << mom[0]    << '\n';
	    std::cout << "TIME= " << time << " YMOM        = "   << mom[1]    << '\n';
	    std::cout << "TIME= " << time << " ZMOM        = "   << mom[2]    << '\n';
	    std::cout << "TIME= " << time << " RHO*e       = "   << rho_e     << '\n';
	    std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
	    if (parent->NumDataLogs() > 0 ) {

	       std::ostream& data_log1 = parent->DataLog(0);

	       if (data_log1.good()) {

		  if (time == 0.0) {
		      data_log1 << std::setw(datwidth) <<  "          time";
		      data_log1 << std::setw(datwidth) <<  "          mass";
		      data_log1 << std::setw(datwidth) <<  "          xmom";
		      data_log1 << std::setw(datwidth) <<  "          ymom";
		      data_log1 << std::setw(datwidth) <<  "          zmom";
		      data_log1 << std::setw(datwidth) <<  "         rho_K";
		      data_log1 << std::setw(datwidth) <<  "         rho_E";
		      data_log1 << std::endl;
		  }

		  // Write the quantities at this time
		  data_log1 << std::setw(datwidth) <<  time;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mass;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[0];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[1];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[2];
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_e;
		  data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_E;
		  data_log1 << std::endl;

	       }

	    }
	    
	    if (show_center_of_mass) {
	        for (int i = 0; i <= 2; i++) {
		  com[i]     = com[i] / mass;
		  com_vel[i] = mom[i] / mass;
		}

		std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com[0]     << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_vel[0] << '\n';

		std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com[1]     << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_vel[1] << '\n';

		std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com[2]     << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_vel[2] << '\n';
	    }
	}
#ifdef BL_LAZY
	});
#endif
    }
}
