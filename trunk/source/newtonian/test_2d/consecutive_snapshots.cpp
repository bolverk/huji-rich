#include "consecutive_snapshots.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots(Trigger* trigger,
					   Index2FileName* i2fn):
  trigger_(trigger), i2fn_(i2fn), counter_(0) {}

void ConsecutiveSnapshots::operator()(hdsim const& sim)
{
  if((*trigger_.get())(sim)){
    write_snapshot_to_hdf5(sim, (*i2fn_.get())(counter_));
    ++counter_;
  }
}
