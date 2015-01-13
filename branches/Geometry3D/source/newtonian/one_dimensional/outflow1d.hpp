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
  
  Conserved CalcFlux(vector<double> const& Vertices,
		     vector<Primitive> const& Cells,
		     RiemannSolver const& rs,
		     vector<double> const& vertex_veclocity,
		     int i) const;
};

#endif // OUTFLOW_1D_HPP
