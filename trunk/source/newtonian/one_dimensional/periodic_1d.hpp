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

  Conserved CalcFlux(vector<double> const& Vertices, 
		     vector<Primitive> const& Cells,
		     RiemannSolver const& rs, 
		     vector<double> const& vertex_velocity,
		     int i) const;
};

#endif // PERIODIC_1D_HPP
