#include "consecutive_snapshots.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots(double dt,
					   Index2FileName* i2fn,
					   double init_time,
					   int counter):
  next_time_(init_time),
  dt_(dt),
  counter_(counter),
  i2fn_(i2fn) {}

void ConsecutiveSnapshots::operator()(hdsim const& sim)
{
  if(sim.GetTime()>next_time_){
    
    write_snapshot_to_hdf5(sim, (*i2fn_.get())(counter_));
 
    next_time_ += dt_;
    ++counter_;
  }
}
