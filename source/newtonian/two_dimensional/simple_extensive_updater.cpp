#include "simple_extensive_updater.hpp"

namespace 
{
  bool bracketed(int low, int arg, int high)
  {
    return arg>=low && high>arg;
  }
}

void SimpleExtensiveUpdater::operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& /*pg*/,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   const vector<ComputationalCell>& /*cells*/,
   vector<Extensive>& extensives) const
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  for(size_t i=0;i<edge_list.size();++i){
    const Edge& edge = edge_list[i];
    const Extensive delta = dt*cd.areas[i]*fluxes[i];
    if(bracketed(0,edge.neighbors.first,tess.GetPointNo()))
      extensives[static_cast<size_t>(edge.neighbors.first)] -=
	delta;
    if(bracketed(0,edge.neighbors.second,tess.GetPointNo()))
      extensives[static_cast<size_t>(edge.neighbors.second)] +=
	delta;
  }
}
