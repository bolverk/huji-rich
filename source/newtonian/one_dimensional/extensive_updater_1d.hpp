#ifndef EXTENSIVE_UPDATER_1D_HPP
#define EXTENSIVE_UPDATER_1D_HPP 1

#include <vector>
#include "../two_dimensional/extensive.hpp"
#include "physical_geometry_1d.hpp"
#include "simulation_state_1d.hpp"

using std::vector;

class ExtensiveUpdater1D
{
public:

  virtual void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   double dt,
   vector<Extensive>& extensives) const = 0;

  virtual ~ExtensiveUpdater1D(void);
};

#endif // EXTENSIVE_UPDATER_1D_HPP
