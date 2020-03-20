#include "simple_extensive_updater_1d.hpp"

SimpleExtensiveUpdater1D::SimpleExtensiveUpdater1D(void) {}

void SimpleExtensiveUpdater1D::operator()
(const vector<Extensive>& fluxes,
 const PhysicalGeometry1D& pg,
 const SimulationState1D& ss,
 double dt,
 vector<Extensive>& extensives) const
{
  for(size_t i=0;i<extensives.size();++i){
    extensives.at(i) += dt*pg.calcArea(ss.getVertices().at(i))*fluxes.at(i);
    extensives.at(i) -= dt*pg.calcArea(ss.getVertices().at(i+1))*fluxes.at(i+1);
  }
}
