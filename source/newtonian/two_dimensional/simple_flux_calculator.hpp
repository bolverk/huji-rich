/*! \file simple_flux_calculator.hpp
  \author Almog Yalinewich
  \brief Simple flux calculator
 */

#ifndef SIMPLE_FLUX_CALCULATOR_HPP
#define SIMPLE_FLUX_CALCULATOR_HPP 1

#include "flux_calculator.hpp"
#include "../common/riemann_solver.hpp"

//! \brief Simple flux calculator
class SimpleFluxCalculator: public FluxCalculator
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
   */
  SimpleFluxCalculator(const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const double time,
   const double dt) const;

private:
  const RiemannSolver& rs_;

  Conserved calcHydroFlux(const Tessellation& tess,
			  const vector<Vector2D>& point_velocities,
			  const vector<ComputationalCell>& cells,
			  const EquationOfState& eos,
			  const size_t i) const;
};

#endif // SIMPLE_FLUX_CALCULATOR_HPP
