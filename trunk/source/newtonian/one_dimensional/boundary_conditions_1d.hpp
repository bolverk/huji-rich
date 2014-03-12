/*! \file boundary_conditions_1d.hpp
  \brief Base class for boundary conditions
  \author Almog Yalinewich
 */

#ifndef BOUNDARY_CONDITIONS_1D_HPP
#define BOUNDARY_CONDITIONS_1D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"
#include "../common/riemann_solver.hpp"

using std::vector;

//! \brief Base class for boundary conditions
class BoundaryConditions1D
{
public:
  /*! \brief Calculates the flux at the boundaries
    \param Vertices Position of the vertices
    \param Cells Hydrodynamic variables
    \param rs Riemann solver
    \param vertex_velocity Velocity of the vertex
    \param i Vertex index
    \return Flux at the boundary
   */
  virtual Conserved CalcFlux
  (vector<double> const& Vertices,
   vector<Primitive> const& Cells,
   RiemannSolver const& rs, 
   vector<double> const& vertex_velocity,
   int i) const= 0;

  virtual ~BoundaryConditions1D(void);
};

#endif // BOUNDARY_CONDITIONS_1D_HPP
