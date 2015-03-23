#ifndef SIMPLE_FLUX_CALCULATOR_HPP
#define SIMPLE_FLUX_CALCULATOR_HPP 1

#include "flux_calculator.hpp"
#include "../common/riemann_solver.hpp"

class SimpleFluxCalculator: public FluxCalculator
{
public:

  SimpleFluxCalculator(const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& face_velocities,
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
