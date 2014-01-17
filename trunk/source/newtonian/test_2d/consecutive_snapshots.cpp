#include "consecutive_snapshots.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots(double dt,double init_time,int counter):
  next_time_(init_time),last_time_(init_time),
  dt_(dt),
  counter_(counter) {}

void ConsecutiveSnapshots::diagnose(hdsim const& sim)
{
  write_number(sim.GetTime(),"time.txt");
  write_number(sim.GetTime()-last_time_,"dt.txt");
  last_time_=sim.GetTime();
  if(sim.GetTime()>next_time_){
    next_time_ += dt_;

    write_snapshot_to_hdf5(sim,"snapshot_"+int2str(counter_)+".h5");
    ++counter_;
  }
}
