/*! \file source_term_1d.hpp
  \brief Abstract class for external forces
  \author Almog Yalinewich
*/

#ifndef EXTERNAL_FORCES_1D_HPP
#define EXTERNAL_FORCES_1D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"
#include "simulation_state_1d.hpp"
#include "../two_dimensional/extensive.hpp"

using std::vector;

//! \brief Abstract class for external forces
class SourceTerm1D
{
public:

  /*! \brief Calculates the change in the extensive conserved variables due to external forces
    \param point Point in which the forces will be evaluated
    \param t Simulation time
    \param dt Time step
    \return Value of the external force
   */
  virtual Extensive operator()
  (const SimulationState1D& state,
   size_t point,
   double t,
   double dt) const = 0;

  virtual ~SourceTerm1D(void);
};

#endif // EXTERNAL_FORCES_1D_HPP
