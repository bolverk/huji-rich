#include "main_loop_1d.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/simple_io.hpp"

void simulation1d::main_loop(hdsim1D& sim,
			     TerminationCondition& term_cond,
			     int time_order,
			     DiagnosticsFunction* diag)
{
  while(term_cond.should_continue(sim)){
    if(1==time_order)
      sim.TimeAdvance();
    else if(2==time_order)
      sim.TimeAdvanceRK(2);
    else
      throw UniversalError("Error in 1d main_loop: unsupported time integration order");

    if(diag)
      diag->diagnose(sim);
  }
}

void simulation1d::main_loop(hdsim1D& sim,
			     double final_time,
			     int max_iter,
			     int time_order,
			     string const& time_log)
{
  while(sim.GetTime()<final_time){
    if(1==time_order)
      sim.TimeAdvance();
    else if(2==time_order)
      sim.TimeAdvanceRK(2);
    else
      throw UniversalError("Error in 1d main_loop: unsupported time integration order");
    
    if(time_log!="")
      write_number(sim.GetTime(),time_log);

    // Infinite loop guard
    if(sim.GetCycle()>max_iter)
      throw UniversalError("Too many time advance iterations");
  }
}
