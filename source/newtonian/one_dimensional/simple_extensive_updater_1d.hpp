/*! \file simple_extensive_updater_1d.hpp
  \author Almog Yalinewich
  \brief A simple extensive updater
 */

#include "extensive_updater_1d.hpp"

//! \brief A simple extensive updater
class SimpleExtensiveUpdater1D: public ExtensiveUpdater1D
{
public:

  //! \brief Class constructor
  SimpleExtensiveUpdater1D(void);

  void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   double dt,
   vector<Extensive>& extensives) const;
};
