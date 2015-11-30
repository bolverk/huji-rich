#include "simple_cell_updater.hpp"
#include "../../misc/lazy_list.hpp"

SimpleCellUpdater::SimpleCellUpdater
(const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence):
  sequence_(sequence) {}

SimpleCellUpdater::~SimpleCellUpdater(void)
{
  for(size_t i=0;i<sequence_.size();++i){
    delete sequence_[i].first;
    delete sequence_[i].second;
  }    
}

SimpleCellUpdater::Condition::~Condition(void) {}

SimpleCellUpdater::Action::~Action(void) {}

namespace {

  ComputationalCell regular_update
  (const EquationOfState& eos,
   vector<Extensive>& extensives,
   const ComputationalCell& old,
   const CacheData& cd,
   const size_t index)
  {
    ComputationalCell res;
    Extensive& extensive = extensives[index];
    const double volume = cd.volumes[index];
    res.density = extensive.mass/volume;
    res.velocity = extensive.momentum/extensive.mass;
    const double energy = extensive.energy/extensive.mass -
      0.5*ScalarProd(res.velocity,res.velocity);
	res.tracers.reserve(extensive.tracers.size());
	for (size_t i = 0; i < extensive.tracers.size(); ++i)
	  res.tracers.insert(pair<string,double>((extensive.tracers.begin() + static_cast<int>(i))->first ,
						 (extensive.tracers.begin() + static_cast<int>(i))->second/extensive.mass));
    res.pressure = eos.de2p
      (res.density,
       energy,
       res.tracers);
    res.stickers = old.stickers;
	string entropy = "Entropy";
	if (old.tracers.find(entropy) != old.tracers.end())
	{
		res.tracers[entropy] = eos.dp2s(res.density, res.pressure);
		extensive.tracers[entropy] = res.tracers[entropy] * extensive.mass;
	}
    return res;
  }

  ComputationalCell update_single
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd,
   const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >& sequence,
   const size_t index)
  {
    for(size_t i=0;i<sequence.size();++i){
      if((*sequence[i].first)
	 (tess,
	  pg,
	  eos,
	  extensives,
	  old,
	  cd,
	  index))
	return (*sequence[i].second)
	 (tess,
	  pg,
	  eos,
	  extensives,
	  old,
	  cd,
	  index);
    }
    return regular_update
      (eos,extensives,old.at(index),cd,index);
  }
}

vector<ComputationalCell> SimpleCellUpdater::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   vector<Extensive>& extensives,
   const vector<ComputationalCell>& old,
   const CacheData& cd) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i)
    res[i] = update_single
      (tess,
       pg,
       eos,
       extensives,
       old,
       cd,
       sequence_,
       i);
  return res;
}

HasSticker::HasSticker
(const string& sticker_name):
  sticker_name_(sticker_name) {}

bool HasSticker::operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& /*eos*/,
   const vector<Extensive>& /*extensives*/,
   const vector<ComputationalCell>& cells,
   const CacheData& /*cd*/,
   const size_t index) const
{
  return safe_retrieve(cells.at(index).stickers,sticker_name_);
}

SkipUpdate::SkipUpdate(void) {}

ComputationalCell SkipUpdate::operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
   const EquationOfState& /*eos*/,
   const vector<Extensive>& /*extensives*/,
   const vector<ComputationalCell>& cells,
   const CacheData& /*cd*/,
   const size_t index) const
{
  return cells[index];
}
