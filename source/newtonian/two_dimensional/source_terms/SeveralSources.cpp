#include "SeveralSources.hpp"

SeveralSources::SeveralSources(vector<SourceTerm*> sources):
  sources_(sources) {}

SeveralSources::~SeveralSources(void) {}

vector<Extensive> SeveralSources::operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const vector<Vector2D>& point_velocities,
   const double t) const
{
  vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i){
    res[i].mass = 0;
    res[i].momentum = Vector2D(0,0);
    res[i].energy = 0;
    for(boost::container::flat_map<std::string,double>::const_iterator it =
	  cells.front().tracers.begin();
	it != cells.front().tracers.end(); ++it)
      res[i].tracers[it->first] = 0;
  }
  for(size_t i=0;i<sources_.size();++i){
    const vector<Extensive> diff = (*sources_[i])
      (tess,pg,cd,cells,fluxes,point_velocities,t);
    for(size_t j=0;j<res.size();++j){
      res[j].mass += diff[j].mass;
      res[j].momentum += diff[j].momentum;
      res[j].energy += diff[j].energy;
      for(boost::container::flat_map<std::string,double>::const_iterator it =
	    diff[j].tracers.begin();
	  it != diff[j].tracers.end(); ++it)
	res[j].tracers[it->first] += it->second;
    }
  }
  return res;
}
