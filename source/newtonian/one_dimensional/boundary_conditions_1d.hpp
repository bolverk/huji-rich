/*! \file boundary_conditions_1d.hpp
  \brief Base class for boundary conditions
  \author Almog Yalinewich
 */

#ifndef BOUNDARY_CONDITIONS_1D_HPP
#define BOUNDARY_CONDITIONS_1D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"
#include "../common/equation_of_state.hpp"
#include "../two_dimensional/extensive.hpp"
#include "simulation_state_1d.hpp"

using std::vector;

//! \brief Base class for boundary conditions
class BoundaryConditions1D
{
public:
  /*! \brief Calculates the flux at the boundaries
    \param ss Computational domain and hydro cells
    \param eos Equation of state
    \param rs Riemann solver
    \param vertex_velocity Velocity of the vertex
    \param side False for left boundary, true for right boundary
    \return Flux at the boundary
   */
  virtual Extensive operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const bool side) const= 0;

  virtual ~BoundaryConditions1D(void);
};

#endif // BOUNDARY_CONDITIONS_1D_HPP
