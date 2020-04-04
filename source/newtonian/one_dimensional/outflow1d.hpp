/*! \file outflow1d.hpp
  \brief Outflow boundary conditions
  \author Almog Yalinewich
*/

#ifndef OUTFLOW_1D_HPP
#define OUTFLOW_1D_HPP 1

#include "boundary_conditions_1d.hpp"

//! \brief Outflow boundary conditions
class Outflow: public BoundaryConditions1D
{
public:
  
  Extensive operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs,
   const vector<double>& vertex_veclocity,
   const size_t i) const;
};

#endif // OUTFLOW_1D_HPP
