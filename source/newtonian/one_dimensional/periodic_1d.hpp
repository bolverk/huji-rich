/*! \file periodic_1d.hpp
  \brief Periodic boundary conditions
  \author Almog Yalinewich
 */

#ifndef PERIODIC_1D_HPP
#define PERIODIC_1D_HPP 1

#include "boundary_conditions_1d.hpp"

/*! \brief Periodic boundary conditions
 */
class Periodic1D: public BoundaryConditions1D
{
public:

  Extensive operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const bool side) const;
};

#endif // PERIODIC_1D_HPP
