#include <cstdio>

#include "AMReX_LevelBld.H"
#include <AMReX_ParmParse.H>
#include "Castro.H"
#include "Castro_F.H"
#include "AMReX_buildInfo.H"

using std::string;
using namespace amrex;

static Box the_same_box (const Box& b) { return b; }

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
  };

static int norm_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD
  };

static int tang_vel_bc[] =
  {
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
  };

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  for (int i = 0; i < BL_SPACEDIM; i++)
    {
      bc.setLo(i,scalar_bc[lo_bc[i]]);
      bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,norm_vel_bc[lo_bc[0]]);
  bc.setHi(0,norm_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,norm_vel_bc[lo_bc[1]]);
  bc.setHi(1,norm_vel_bc[hi_bc[1]]);
  bc.setLo(2,tang_vel_bc[lo_bc[2]]);
  bc.setHi(2,tang_vel_bc[hi_bc[2]]);
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();
  bc.setLo(0,tang_vel_bc[lo_bc[0]]);
  bc.setHi(0,tang_vel_bc[hi_bc[0]]);
  bc.setLo(1,tang_vel_bc[lo_bc[1]]);
  bc.setHi(1,tang_vel_bc[hi_bc[1]]);
  bc.setLo(2,norm_vel_bc[lo_bc[2]]);
  bc.setHi(2,norm_vel_bc[hi_bc[2]]);
}

void
Castro::variableSetUp ()
{

  // Castro::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job


  // initialize the start time for our CPU-time tracker
  startCPUTime = ParallelDescriptor::second();


  // Output the git commit hashes used to build the executable.

  if (ParallelDescriptor::IOProcessor()) {

    const char* castro_hash  = buildInfoGetGitHash(1);
    const char* amrex_hash  = buildInfoGetGitHash(2);

    if (strlen(castro_hash) > 0) {
      std::cout << "\n" << "Castro git describe: " << castro_hash << "\n";
    }
    if (strlen(amrex_hash) > 0) {
      std::cout << "AMReX git describe: " << amrex_hash << "\n";
    }

    std::cout << "\n";
  }

  BL_ASSERT(desc_lst.size() == 0);

  // Get options, set phys_bc
  read_params();

  // Initialize the runtime parameters for any of the external
  // microphysics
  extern_init();

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

#ifdef NUM_ADV
  NumAdv = NUM_ADV;
#else
  NumAdv = 0;
#endif

  if (NumAdv > 0)
    {
      FirstAdv = cnt;
      cnt += NumAdv;
    }

  int dm = BL_SPACEDIM;

  // Get the number of species from the network model.
  ca_get_num_spec(&NumSpec);

  if (NumSpec > 0)
    {
      FirstSpec = cnt;
      cnt += NumSpec;
    }

  // Get the number of auxiliary quantities from the network model.
  ca_get_num_aux(&NumAux);

  if (NumAux > 0)
    {
      FirstAux = cnt;
      cnt += NumAux;
    }

  NUM_STATE = cnt;

  // Define NUM_GROW from the F90 module.
  ca_get_method_params(&NUM_GROW);

  const Real run_strt = ParallelDescriptor::second() ;

  // Read in the input values to Fortran.

  ca_set_castro_method_params();

  ca_set_method_params(dm, Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec, FirstAux,
		       NumAdv);

  // Get the number of primitive variables from Fortran.

  ca_get_qvar(&QVAR);
  ca_get_nqaux(&NQAUX);

  Real run_stop = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_stop,ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor())
    std::cout << "\nTime in ca_set_method_params: " << run_stop << '\n' ;

  int coord_type = Geometry::Coord();

  // Get the center variable from the inputs and pass it directly to Fortran.
  Array<Real> center(BL_SPACEDIM, 0.0);
  ParmParse ppc("castro");
  ppc.queryarr("center",center,0,BL_SPACEDIM);

  ca_set_problem_params(dm,phys_bc.lo(),phys_bc.hi(),
			Interior,Inflow,Outflow,Symmetry,SlipWall,NoSlipWall,coord_type,
			Geometry::ProbLo(),Geometry::ProbHi(),center.dataPtr());

  // Read in the parameters for the tagging criteria
  // and store them in the Fortran module.

  int probin_file_length = probin_file.length();
  Array<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  ca_get_tagging_params(probin_file_name.dataPtr(),&probin_file_length);

  Interpolater* interp = &cell_cons_interp;

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = 0;

  store_in_checkpoint = true;
  desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,ngrow_state,NUM_STATE,
			 interp,state_data_extrap,store_in_checkpoint);

  Array<BCRec>       bcs(NUM_STATE);
  Array<std::string> name(NUM_STATE);

  BCRec bc;
  cnt = 0;
  set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
  cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
  cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
  cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
  cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
  cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
  cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";

  for (int i=0; i<NumAdv; ++i)
    {
      char buf[64];
      sprintf(buf, "adv_%d", i);
      cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = string(buf);
    }

  // Get the species names from the network model.
  std::vector<std::string> spec_names;
  for (int i = 0; i < NumSpec; i++) {
    int len = 20;
    Array<int> int_spec_names(len);
    // This call return the actual length of each string in "len"
    ca_get_spec_names(int_spec_names.dataPtr(),&i,&len);
    char char_spec_names[len+1];
    for (int j = 0; j < len; j++)
      char_spec_names[j] = int_spec_names[j];
    char_spec_names[len] = '\0';
    spec_names.push_back(std::string(char_spec_names));
  }

  if ( ParallelDescriptor::IOProcessor())
    {
      std::cout << NumSpec << " Species: " << std::endl;
      for (int i = 0; i < NumSpec; i++)
	std::cout << spec_names[i] << ' ' << ' ';
      std::cout << std::endl;
    }

  for (int i=0; i<NumSpec; ++i)
    {
      cnt++;
      set_scalar_bc(bc,phys_bc);
      bcs[cnt] = bc;
      name[cnt] = "rho_" + spec_names[i];
    }

  // Get the auxiliary names from the network model.
  std::vector<std::string> aux_names;
  for (int i = 0; i < NumAux; i++) {
    int len = 20;
    Array<int> int_aux_names(len);
    // This call return the actual length of each string in "len"
    ca_get_aux_names(int_aux_names.dataPtr(),&i,&len);
    char char_aux_names[len+1];
    for (int j = 0; j < len; j++)
      char_aux_names[j] = int_aux_names[j];
    char_aux_names[len] = '\0';
    aux_names.push_back(std::string(char_aux_names));
  }

  if ( ParallelDescriptor::IOProcessor())
    {
      std::cout << NumAux << " Auxiliary Variables: " << std::endl;
      for (int i = 0; i < NumAux; i++)
	std::cout << aux_names[i] << ' ' << ' ';
      std::cout << std::endl;
    }

  for (int i=0; i<NumAux; ++i)
    {
      cnt++;
      set_scalar_bc(bc,phys_bc);
      bcs[cnt] = bc;
      name[cnt] = "rho_" + aux_names[i];
    }

  desc_lst.setComponent(State_Type,
			Density,
			name,
			bcs,
			BndryFunc(ca_denfill,ca_hypfill));

  num_state_type = desc_lst.size();

  //
  // DEFINE DERIVED QUANTITIES
  //

  //
  // Pressure
  //
  derive_lst.add("pressure",IndexType::TheCellType(),1,ca_derpres,the_same_box);
  derive_lst.addComponent("pressure",desc_lst,State_Type,Density,NUM_STATE);

  //
  // Velocities
  //
  derive_lst.add("x_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
  derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("x_velocity",desc_lst,State_Type,Xmom,1);

  derive_lst.add("y_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
  derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("y_velocity",desc_lst,State_Type,Ymom,1);

  derive_lst.add("z_velocity",IndexType::TheCellType(),1,ca_dervel,the_same_box);
  derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,1);
  derive_lst.addComponent("z_velocity",desc_lst,State_Type,Zmom,1);

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
