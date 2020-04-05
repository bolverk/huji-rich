/*! \file simple_flux_calculator_1d.hpp
  \author Almog Yalinewich
  \brief Simple flux calculator
 */

#ifndef SIMPLE_FLUX_CALCULATOR_1D_HPP
#define SIMPLE_FLUX_CALCULATOR_1D_HPP 1

#include "flux_calculator_1d.hpp"
#include "../common/riemann_solver.hpp"
#include "spatial_reconstruction1d.hpp"
#include "boundary_conditions_1d.hpp"

//! \brief Simple flux calculator
class SimpleFluxCalculator1D: public FluxCalculator1D
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
    \param interp Interpolation scheme
    \param bc Boundary conditions
   */
  SimpleFluxCalculator1D
  (const RiemannSolver& rs,
   const SpatialReconstruction1D& interp,
   const BoundaryConditions1D& bc);

  vector<Extensive> operator()
  (const SimulationState1D& ss,
   const vector<double>& vertex_velocity,
   const EquationOfState& eos,
   const double dt) const;

private:
  const RiemannSolver& rs_;
  const SpatialReconstruction1D& interp_;
  const BoundaryConditions1D& bc_;
};

#endif // SIMPLE_FLUX_CALCULATOR_1D_HPP
