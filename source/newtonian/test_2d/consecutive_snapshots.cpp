#include "consecutive_snapshots.hpp"
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include "../two_dimensional/hdf5_diagnostics.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots(Trigger* trigger,
					   Index2FileName* i2fn):
  trigger_(trigger), i2fn_(i2fn), counter_(0), appendices_() {}

ConsecutiveSnapshots::ConsecutiveSnapshots(Trigger* trigger,
					   Index2FileName* i2fn,
					   const vector<DiagnosticAppendix*>& appendices):
  trigger_(trigger), i2fn_(i2fn), counter_(0), appendices_(appendices) {}

void ConsecutiveSnapshots::operator()(hdsim const& sim)
{
  if((*trigger_.get())(sim)){
    write_snapshot_to_hdf5(sim, (*i2fn_.get())(counter_), appendices_);
    ++counter_;
  }
}

ConsecutiveSnapshots::~ConsecutiveSnapshots(void)
{
  for(size_t i=0;i<appendices_.size();++i)
    delete appendices_[i];
}
