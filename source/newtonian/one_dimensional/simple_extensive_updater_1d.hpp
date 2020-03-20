#include "extensive_updater_1d.hpp"

class SimpleExtensiveUpdater1D: public ExtensiveUpdater1D
{
public:

  SimpleExtensiveUpdater1D(void);

  void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   double dt,
   vector<Extensive>& extensives) const;
};
