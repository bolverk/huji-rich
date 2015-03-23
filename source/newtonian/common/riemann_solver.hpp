/*! \file riemann_solver.hpp
  \brief Base class for Riemann solver
  \author Almog Yalinewich
 */

#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP 1

#include "hydrodynamic_variables.hpp"

//! \brief Base class for Riemann solver
class RiemannSolver
{
public:

  /*! \brief Solve Riemann porblme
    \param left Primitive variables on the left side
    \param right Primitive variables on the right side
    \param velocity Velocity of the vertex
    \return Corrected fluxes
   */
  virtual Conserved operator()
  (Primitive const& left, Primitive const& right, double velocity) const = 0;

  virtual ~RiemannSolver(void);
};

#endif
