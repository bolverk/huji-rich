#include "consecutive_snapshots.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots(double dt):
  next_time_(0),
  dt_(dt),
  counter_(0) {}

void ConsecutiveSnapshots::diagnose(hdsim const& sim)
{
  write_number(sim.GetTime(),"time.txt");

  if(sim.GetTime()>next_time_){
    next_time_ += dt_;

    write_snapshot_to_hdf5(sim,"snapshot_"+int2str(counter_)+".h5");
    ++counter_;
  }
}
