/*! \file rigid_wall_1d.hpp
  \brief Rigid wall boundary conditions
  \author Almog Yalinewich
 */

#ifndef RIGID_WALL_1D_HPP
#define RIGID_WALL_1D_HPP 1

#include "boundary_conditions_1d.hpp"

//! \brief Rigid wall boundary condition in 1d
class RigidWall1D: public BoundaryConditions1D
{
public:

  Extensive operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const size_t i) const;
};

#endif // RIGID_WALL_1D_HPP
