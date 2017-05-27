#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::apply_source_to_state(MultiFab& state, MultiFab& source, Real dt)
{

  MultiFab::Saxpy(state, dt, source, 0, 0, NUM_STATE, 0);

}

void
Castro::time_center_source_terms(MultiFab& S_new, MultiFab& src_old, MultiFab &src_new, Real dt)
{
  BL_PROFILE("Castro::time_center_source_terms()");

  // Subtract off half of the old source term, and add half of the new.

  MultiFab::Saxpy(S_new,-0.5*dt,src_old,0,0,S_new.nComp(),0);
  MultiFab::Saxpy(S_new, 0.5*dt,src_new,0,0,S_new.nComp(),0);
}

bool
Castro::source_flag(int src)
{
    switch(src) {

#ifdef SPONGE
    case sponge_src:
	if (do_sponge)
	    return true;
	else
	    return false;
#endif

    case ext_src:
	if (add_ext_src)
	    return true;
	else
	    return false;

#ifdef GRAVITY
    case grav_src:
	if (do_grav)
	    return true;
	else
	    return false;
#endif

    default:
	return false;

    } // end switch
}

void
Castro::do_old_sources(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

    // Construct the old-time sources.

    for (int n = 0; n < num_src; ++n)
	construct_old_source(n, time, dt, amr_iteration, amr_ncycle,
			     sub_iteration, sub_ncycle);


    if (do_ctu) {
      // Apply the old-time sources directly to the new-time state,
      // S_new -- note that this addition is for full dt, since we
      // will do a predictor-corrector on the sources to allow for
      // state-dependent sources.

      // note: this is not needed for MOL, since the sources will
      // be applied as part of the overall integration

      MultiFab& S_new = get_new_data(State_Type);

      for (int n = 0; n < num_src; ++n)
	if (source_flag(n))
	  apply_source_to_state(S_new, *old_sources[n], dt);
    }

    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (print_update_diagnostics) {
      bool is_new = false;
      print_all_source_changes(dt, is_new);
    }

}

void
Castro::do_new_sources(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

    MultiFab& S_new = get_new_data(State_Type);

    // For the new-time source terms, we have an option for how to proceed.
    // We can either construct all of the old-time sources using the same
    // state that comes out of the hydro update, or we can evaluate the sources
    // one by one and apply them as we go.

    if (update_state_between_sources) {

	for (int n = 0; n < num_src; ++n) {
	    construct_new_source(n, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);
	    if (source_flag(n)) {
		apply_source_to_state(S_new, *new_sources[n], dt);
		clean_state(S_new);
	    }
	}

    } else {

	// Construct the new-time source terms.

	for (int n = 0; n < num_src; ++n)
	    construct_new_source(n, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

	// Apply the new-time sources to the state.

	for (int n = 0; n < num_src; ++n)
	    if (source_flag(n))
		apply_source_to_state(S_new, *new_sources[n], dt);

	clean_state(S_new);

    }

    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (print_update_diagnostics) {
      bool is_new = true;
      print_all_source_changes(dt, is_new);
    }


}

void
Castro::construct_old_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
	construct_old_sponge_source(time, dt);
	break;
#endif

    case ext_src:
	construct_old_ext_source(time, dt);
	break;

#ifdef GRAVITY
    case grav_src:
	construct_old_gravity_source(time, dt);
	break;
#endif

    default:
	break;

    } // end switch
}

void
Castro::construct_new_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
	construct_new_sponge_source(time, dt);
	break;
#endif

    case ext_src:
	construct_new_ext_source(time, dt);
	break;

#ifdef GRAVITY
    case grav_src:
	construct_new_gravity_source(time, dt);
	break;
#endif

    default:
	break;

    } // end switch
}

// Evaluate diagnostics quantities describing the effect of an
// update on the state. The optional parameter local determines
// whether we want to do this calculation globally over all processes
// or locally just on this processor. The latter is useful if you
// are evaluating the contribution from multiple source changes at once
// and want to amortize the cost of the parallel reductions.
// Note that the resultant output is volume-weighted.

Array<Real>
Castro::evaluate_source_change(MultiFab& source, Real dt, bool local)
{

  Array<Real> update(source.nComp(), 0.0);

  // Create a temporary array which will hold a single component
  // at a time of the volume-weighted source.

  MultiFab weighted_source(source.boxArray(), source.DistributionMap(), 1, 0);

  for (int n = 0; n < source.nComp(); ++n) {

    weighted_source.setVal(0.0);

    // Fill weighted_source with source x volume.

    MultiFab::AddProduct(weighted_source, source, n, volume, 0, 0, 1, 0);

    update[n] = weighted_source.sum(0, local) * dt;

  }

  return update;

}

// Print the change due to a given source term update.
// We assume here that the input array is lined up with
// the NUM_STATE components of State_Type because we are
// interested in printing changes to energy, mass, etc.

void
Castro::print_source_change(Array<Real> update)
{

  BL_ASSERT(update.size() == NUM_STATE);

  if (ParallelDescriptor::IOProcessor()) {

    std::cout << "       mass added: " << update[Density] << std::endl;
    std::cout << "       xmom added: " << update[Xmom] << std::endl;
#if (BL_SPACEDIM >= 2)
    std::cout << "       ymom added: " << update[Ymom] << std::endl;
#endif
#if (BL_SPACEDIM == 3)
    std::cout << "       zmom added: " << update[Zmom] << std::endl;
#endif
    std::cout << "       eint added: " << update[Eint] << std::endl;
    std::cout << "       ener added: " << update[Eden] << std::endl;

    std::cout << std::endl;

  }


}

// For the old-time or new-time sources update, evaluate the change in the state
// for all source terms, then pring the results.

void
Castro::print_all_source_changes(Real dt, bool is_new)
{

  Array< Array<Real> > summed_updates;

  summed_updates.resize(num_src);

  bool local = true;

  for (int n = 0; n < num_src; ++n) {

    if (!source_flag(n)) continue;

    MultiFab& source = is_new ? *new_sources[n] : *old_sources[n];

    summed_updates[n] = evaluate_source_change(source, dt, local);

  }

#ifdef BL_LAZY
  Lazy::QueueReduction( [=] () mutable {
#endif
      if (coalesce_update_diagnostics) {

	  Array<Real> coalesced_update(NUM_STATE, 0.0);

	  for (int n = 0; n < num_src; ++n) {
	      if (!source_flag(n)) continue;

	      for (int s = 0; s < NUM_STATE; ++s) {
		  coalesced_update[s] += summed_updates[n][s];
	      }
	  }

	  ParallelDescriptor::ReduceRealSum(coalesced_update.dataPtr(), NUM_STATE, ParallelDescriptor::IOProcessorNumber());

	  std::string time = is_new ? "new" : "old";

	  if (ParallelDescriptor::IOProcessor())
	      std::cout << std::endl << "  Contributions to the state from the " << time << "-time sources:" << std::endl;

	  print_source_change(coalesced_update);

      } else {

	  for (int n = 0; n < num_src; ++n) {

	      if (!source_flag(n)) continue;

	      ParallelDescriptor::ReduceRealSum(summed_updates[n].dataPtr(), NUM_STATE, ParallelDescriptor::IOProcessorNumber());

	      std::string time = is_new ? "new" : "old";

	      if (ParallelDescriptor::IOProcessor())
		  std::cout << std::endl << "  Contributions to the state from the " << time << "-time " << source_names[n] << " source:" << std::endl;

	      print_source_change(summed_updates[n]);

	  }

      }

#ifdef BL_LAZY
    });
#endif

} 

// Obtain the sum of all source terms.

void
Castro::sum_of_sources(MultiFab& source)
{

  // this computes advective_source + 1/2 (old source + new source)
  // 
  // Note: the advective source is defined as -div{F}
  //
  // the time-centering is accomplished since new source is defined
  // to be 1/2 (new source - old source) generally.

  int ng = source.nGrow();

  source.setVal(0.0);

  for (int n = 0; n < num_src; ++n)
      MultiFab::Add(source, *old_sources[n], 0, 0, NUM_STATE, ng);

  MultiFab::Add(source, hydro_source, 0, 0, NUM_STATE, ng);

  for (int n = 0; n < num_src; ++n)
      MultiFab::Add(source, *new_sources[n], 0, 0, NUM_STATE, ng);

}

// Obtain the effective source term due to reactions on the primitive variables.

