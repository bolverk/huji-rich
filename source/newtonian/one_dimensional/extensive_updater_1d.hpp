/*! \file extensive_updater_1d.hpp
  \author Almog Yalinewich
  \brief Extensive updater abstract class
 */

#ifndef EXTENSIVE_UPDATER_1D_HPP
#define EXTENSIVE_UPDATER_1D_HPP 1

#include <vector>
#include "../two_dimensional/extensive.hpp"
#include "physical_geometry_1d.hpp"
#include "simulation_state_1d.hpp"

using std::vector;

//! \brief Method for updating the extensive variables
class ExtensiveUpdater1D
{
public:

  /*! \brief Updates the extensive variabels
    \param fluxes Hydrodynamic fluxes
    \param pg Physical geometry
    \param ss Computational domains and hydro cells
    \param dt Time step
    \param extensives Extensive variables to be updated
   */
  virtual void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   const double dt,
   vector<Extensive>& extensives) const = 0;

  virtual ~ExtensiveUpdater1D(void);
};

#endif // EXTENSIVE_UPDATER_1D_HPP
