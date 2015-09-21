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
   const vector<Extensive>& extensives,
   const ComputationalCell& old,
   const CacheData& cd,
   const size_t index)
  {
    ComputationalCell res;
    const Extensive& extensive = extensives[index];
    const double volume = cd.volumes[index];
    res.density = extensive.mass/volume;
    res.velocity = extensive.momentum/extensive.mass;
    const double energy = extensive.energy/extensive.mass -
      0.5*ScalarProd(res.velocity,res.velocity);
    for(boost::container::flat_map<string,double>::const_iterator it=
	  extensive.tracers.begin();
	it!=extensive.tracers.end();++it)
      res.tracers[it->first] = it->second/extensive.mass;
    res.pressure = eos.de2p
      (res.density,
       energy,
       res.tracers);
    res.stickers = old.stickers;
    return res;
  }

  ComputationalCell update_single
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
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
   const vector<Extensive>& extensives,
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
    /*
    const double volume = cd.volumes[i];
    res[i].density = extensives[i].mass/volume;
    res[i].velocity = extensives[i].momentum / extensives[i].mass;
    const double energy = extensives[i].energy/extensives[i].mass - 
      0.5*ScalarProd(res[i].velocity, res[i].velocity);
    for(boost::container::flat_map<std::string,double>::const_iterator it =
	  extensives[i].tracers.begin();
	it!=extensives[i].tracers.end();++it)
      res[i].tracers[it->first] = it->second/extensives[i].mass;
    res[i].pressure = eos.de2p(res[i].density, energy, res[i].tracers);
    */
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
  return cells.at(index).stickers.find(sticker_name_)->second;
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
