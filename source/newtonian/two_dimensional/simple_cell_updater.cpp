#include "simple_cell_updater.hpp"
#include "../../misc/lazy_list.hpp"

SimpleCellUpdater::SimpleCellUpdater(void) {}

vector<ComputationalCell> SimpleCellUpdater::operator()
  (const Tessellation& /*tess*/,
   const PhysicalGeometry& /*pg*/,
 const EquationOfState& eos,
 const vector<Extensive>& extensives,
 const vector<ComputationalCell>& old,
 const CacheData& cd) const
{
  vector<ComputationalCell> res = old;
  for(size_t i=0;i<extensives.size();++i){
    const double volume = cd.volumes[i];
    res[i].density = extensives[i].mass/volume;
    res[i].velocity = extensives[i].momentum / extensives[i].mass;
    const double energy = extensives[i].energy/extensives[i].mass - 
      0.5*ScalarProd(res[i].velocity, res[i].velocity);
    for(std::map<std::string,double>::const_iterator it =
	  extensives[i].tracers.begin();
	it!=extensives[i].tracers.end();++it)
      res[i].tracers[it->first] = it->second/extensives[i].mass;
    res[i].pressure = eos.de2p(res[i].density, energy, res[i].tracers);
  }
  return res;
}
