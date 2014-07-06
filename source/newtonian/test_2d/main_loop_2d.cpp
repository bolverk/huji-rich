#include "main_loop_2d.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/simple_io.hpp"
#include "../../newtonian/two_dimensional/diagnostics.hpp"

DiagnosticFunction::~DiagnosticFunction(void) {}

WriteTime::WriteTime(string const& fname):
fname_(fname) {}

void WriteTime::operator()(hdsim const& sim)
{
  write_number(sim.GetTime(),fname_);
}

TerminationCondition::~TerminationCondition(void) {}

SafeTimeTermination::SafeTimeTermination
	(double termination_time,
	int max_cycles):
termination_time_(termination_time),
	max_cycles_(max_cycles) {}

bool SafeTimeTermination::operator()(hdsim const& sim)
{
	if(sim.GetCycle()>max_cycles_)
		throw UniversalError("Error in SafeTimeTermination: too many iterations");

	return sim.GetTime()<termination_time_;
}

CycleTermination::CycleTermination(int max_cycles):
max_cycles_(max_cycles) {}

bool CycleTermination::operator()(hdsim const& sim)
{
	return sim.GetCycle()<max_cycles_;
}

void simulation2d::main_loop(hdsim& sim,
	TerminationCondition& term_cond,
	int time_order,
	DiagnosticFunction* diagfunc)
{
	while(term_cond(sim))
	{
		try
		{
			if(1==time_order)
				sim.TimeAdvance();
			else if(2==time_order)
				sim.TimeAdvance2Mid();
			else
				throw UniversalError("Error in 2d main loop: unsupported time integration order");
		}
		catch(UniversalError const& eo)
		{
		  DisplayError(eo);
		}
		if(diagfunc)
		  (*diagfunc)(sim);
	}
}
