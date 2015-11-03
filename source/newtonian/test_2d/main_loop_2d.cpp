#include "main_loop_2d.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/simple_io.hpp"
#include "../../newtonian/two_dimensional/diagnostics.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

DiagnosticFunction::~DiagnosticFunction(void) {}

WriteTime::WriteTime(string const& fname):
  fname_(fname) {}

void WriteTime::operator()(hdsim const& sim)
{
  write_number(sim.getTime(),fname_);
}

TerminationCondition::~TerminationCondition(void) {}

SafeTimeTermination::SafeTimeTermination
(double termination_time,
 int max_cycles):
  termination_time_(termination_time),
  max_cycles_(max_cycles) {}

bool SafeTimeTermination::operator()(hdsim const& sim)
{
  if(sim.getCycle()>max_cycles_)
    throw UniversalError("Error in SafeTimeTermination: too many iterations");

  return sim.getTime()<termination_time_;
}

CycleTermination::CycleTermination(int max_cycles):
  max_cycles_(max_cycles) {}

bool CycleTermination::operator()(hdsim const& sim)
{
  return sim.getCycle()<max_cycles_;
}

Manipulate::~Manipulate(void) {}

void simulation2d::main_loop(hdsim& sim,
			     TerminationCondition& term_cond,
			     void (hdsim::*time_advance_method)(void),
			     DiagnosticFunction* diagfunc,
			     Manipulate* manipulate)
{
  while(term_cond(sim))
    {
      try
	{
	  (sim.*time_advance_method)();
	}
      catch(UniversalError const& eo)
	{
	  DisplayError(eo);
	  throw;
	}
      if(manipulate)
	(*manipulate)(sim);
      if(diagfunc)
	(*diagfunc)(sim);
    }
}

WriteData::WriteData(string const& fname) :
	fname_(fname) {}

void WriteData::operator()(hdsim const& sim)
{
	write_snapshot_to_hdf5(sim, fname_);
}
