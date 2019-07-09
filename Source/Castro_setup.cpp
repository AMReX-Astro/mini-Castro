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
    
void ca_hypfill(BL_FORT_FAB_ARG(state),
                const int dlo[], const int dhi[],
                const amrex::Real dx[], const amrex::Real glo[],
                const amrex::Real* time, const int bc[]) {}

void ca_denfill(BL_FORT_FAB_ARG(state),
                const int dlo[], const int dhi[],
                const amrex::Real dx[], const amrex::Real glo[],
                const amrex::Real* time, const int bc[]) {}

typedef StateDescriptor::BndryFunc BndryFunc;

static
void
set_scalar_bc (BCRec& bc)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, INT_DIR);
        bc.setHi(i, INT_DIR);
    }
}

static
void
set_x_vel_bc(BCRec& bc)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, INT_DIR);
        bc.setHi(i, INT_DIR);
    }
}

static
void
set_y_vel_bc(BCRec& bc)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, INT_DIR);
        bc.setHi(i, INT_DIR);
    }
}

static
void
set_z_vel_bc(BCRec& bc)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, INT_DIR);
        bc.setHi(i, INT_DIR);
    }
}

void
Castro::variableSetUp ()
{

  // Castro::variableSetUp is called in the constructor of Amr.cpp, so
  // it will get called when we start the job.

  BL_ASSERT(desc_lst.size() == 0);

  // Get options, set phys_bc
  read_params();

  // Initialize the network
  network_init();

  //
  // Set number of state variables and pointers to components
  //

  int cnt = 0;
  Density = cnt++;
  Xmom = cnt++;
  Ymom = cnt++;
  Zmom = cnt++;
  Eden = cnt++;
  Eint = cnt++;
  Temp = cnt++;

  NumAdv = 0;

  // Get the number of species from the network model.
  ca_get_num_spec(&NumSpec);

  if (NumSpec > 0) {
      FirstSpec = cnt;
      cnt += NumSpec;
  }

  NumAux = 0;

  NUM_STATE = cnt;

  // Define NUM_GROW from the F90 module.
  ca_get_method_params(&NUM_GROW);

  // Read in the input values to Fortran.

  ca_set_method_params(Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux, NumAdv);

  // Get the number of primitive variables from Fortran.

  ca_get_qvar(&QVAR);
  ca_get_nqaux(&NQAUX);
  ca_get_ngdnv(&NGDNV);

  Interpolater* interp = &cell_cons_interp;

  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = 0;

  store_in_checkpoint = true;
  desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,ngrow_state,NUM_STATE,
			 interp,state_data_extrap,store_in_checkpoint);

//  desc_lst.setDeviceCopy(State_Type, true);

  Vector<BCRec>       bcs(NUM_STATE);
  Vector<std::string> name(NUM_STATE);

  BCRec bc;
  cnt=0; set_scalar_bc(bc); bcs[cnt] = bc; name[cnt] = "density";
  cnt++; set_x_vel_bc(bc);  bcs[cnt] = bc; name[cnt] = "xmom";
  cnt++; set_y_vel_bc(bc);  bcs[cnt] = bc; name[cnt] = "ymom";
  cnt++; set_z_vel_bc(bc);  bcs[cnt] = bc; name[cnt] = "zmom";
  cnt++; set_scalar_bc(bc); bcs[cnt] = bc; name[cnt] = "rho_E";
  cnt++; set_scalar_bc(bc); bcs[cnt] = bc; name[cnt] = "rho_e";
  cnt++; set_scalar_bc(bc); bcs[cnt] = bc; name[cnt] = "Temp";

  // Get the species names from the network model.
  std::vector<std::string> spec_names;
  for (int i = 0; i < NumSpec; i++) {
    int len = 20;
    Vector<int> int_spec_names(len);
    // This call return the actual length of each string in "len"
    ca_get_spec_names(int_spec_names.dataPtr(),&i,&len);
    char char_spec_names[len+1];
    for (int j = 0; j < len; j++)
      char_spec_names[j] = int_spec_names[j];
    char_spec_names[len] = '\0';
    spec_names.push_back(std::string(char_spec_names));
  }

  for (int i=0; i<NumSpec; ++i) {
      cnt++;
      set_scalar_bc(bc);
      bcs[cnt] = bc;
      name[cnt] = "rho_" + spec_names[i];
  }

  desc_lst.setComponent(State_Type,
			Density,
			name,
			bcs,
			BndryFunc(ca_denfill,ca_hypfill));

  num_state_type = desc_lst.size();

  //
  // DEFINE ERROR ESTIMATION QUANTITIES
  //
  ErrorSetUp();

  // method of lines Butcher tableau
#define SECONDORDER_TVD

#ifdef THIRDORDER
  MOL_STAGES = 3;

  a_mol.resize(MOL_STAGES);
  for (int n = 0; n < MOL_STAGES; ++n)
    a_mol[n].resize(MOL_STAGES);

  a_mol[0] = {0,   0, 0};
  a_mol[1] = {0.5, 0, 0};
  a_mol[2] = {-1,  2, 0};

  b_mol = {1./6., 2./3., 1./6.};

  c_mol = {0.0, 0.5, 1};
#endif

#ifdef THIRDORDER_TVD
  MOL_STAGES = 3;

  a_mol.resize(MOL_STAGES);
  for (int n = 0; n < MOL_STAGES; ++n)
    a_mol[n].resize(MOL_STAGES);

  a_mol[0] = {0.0,  0.0,  0.0};
  a_mol[1] = {1.0,  0.0,  0.0};
  a_mol[2] = {0.25, 0.25, 0.0};

  b_mol = {1./6., 1./6., 2./3.};

  c_mol = {0.0, 1.0, 0.5};
#endif

#ifdef SECONDORDER
  MOL_STAGES = 2;

  a_mol.resize(MOL_STAGES);
  for (int n = 0; n < MOL_STAGES; ++n)
    a_mol[n].resize(MOL_STAGES);

  a_mol[0] = {0,   0,};
  a_mol[1] = {0.5, 0,};

  b_mol = {0.0, 1.0};

  c_mol = {0.0, 0.5};
#endif

#ifdef SECONDORDER_TVD
  MOL_STAGES = 2;

  a_mol.resize(MOL_STAGES);
  for (int n = 0; n < MOL_STAGES; ++n)
    a_mol[n].resize(MOL_STAGES);

  a_mol[0] = {0,   0,};
  a_mol[1] = {1.0, 0,};

  b_mol = {0.5, 0.5};

  c_mol = {0.0, 1.0};
#endif

#ifdef FIRSTORDER
  MOL_STAGES = 1;

  a_mol.resize(MOL_STAGES);
  for (int n = 0; n < MOL_STAGES; ++n)
    a_mol[n].resize(MOL_STAGES);

  a_mol[0] = {1};
  b_mol = {1.0};
  c_mol = {0.0};
#endif

}
