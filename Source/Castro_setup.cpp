#include <cstdio>

#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>
#include <Castro.H>
#include <Castro_F.H>

using std::string;
using namespace amrex;

// Hypfill and denfill are empty because we are requiring that
// the boundary conditions are periodic for mini-Castro. However
// the function signatures need to exist so that we can compile.

void hypfill(BL_FORT_FAB_ARG(state),
             const int dlo[], const int dhi[],
             const amrex::Real dx[], const amrex::Real glo[],
             const amrex::Real* time, const int bc[]) {}

void denfill(BL_FORT_FAB_ARG(state),
             const int dlo[], const int dhi[],
             const amrex::Real dx[], const amrex::Real glo[],
             const amrex::Real* time, const int bc[]) {}

typedef StateDescriptor::BndryFunc BndryFunc;

void
Castro::variableSetUp()
{
    BL_PROFILE("Castro::variableSetUp()");

    // Castro::variableSetUp is called in the constructor of Amr.cpp, so
    // it will get called when we start the job.

    BL_ASSERT(desc_lst.size() == 0);

    int bndry_func_thread_safe = 1;
    StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

    // Initialize the network and EOS
    network_init();
    eos_init();

    Interpolater* interp = &cell_cons_interp;

    bool state_data_extrap = false;
    bool store_in_checkpoint;

    int ngrow_state = 0;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,ngrow_state,NUM_STATE,
                           interp,state_data_extrap,store_in_checkpoint);

    Vector<BCRec>       bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);

    BCRec bc;
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, INT_DIR);
        bc.setHi(i, INT_DIR);
    }

    int cnt;
    cnt=0; bcs[cnt] = bc; name[cnt] = "density";
    cnt++; bcs[cnt] = bc; name[cnt] = "xmom";
    cnt++; bcs[cnt] = bc; name[cnt] = "ymom";
    cnt++; bcs[cnt] = bc; name[cnt] = "zmom";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_e";
    cnt++; bcs[cnt] = bc; name[cnt] = "Temp";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_he4";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_c12";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_o16";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_ne20";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_mg24";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_si28";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_s32";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_ar36";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_ca40";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_ti44";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_cr48";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_fe52";
    cnt++; bcs[cnt] = bc; name[cnt] = "rho_ni56";

    desc_lst.setComponent(State_Type, Density, name, bcs, BndryFunc(denfill, hypfill));

    // Update the diagnostic interval.
    ParmParse pp;
    pp.query("diagnostic_interval", diagnostic_interval);
}
