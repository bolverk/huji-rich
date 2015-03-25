/*! \file simple_flux_calculator.hpp
  \author Almog Yalinewich
  \brief Simple flux calculator
 */

#ifndef SIMPLE_FLUX_CALCULATOR_HPP
#define SIMPLE_FLUX_CALCULATOR_HPP 1

#include "flux_calculator_2d.hpp"
#include "../common/riemann_solver.hpp"

Primitive convert_to_primitive(const ComputationalCell& cell,
			       const EquationOfState& eos);

Primitive reflect(const Primitive& p,
		  const Vector2D& axis);

Vector2D remove_parallel_component(const Vector2D& v,
				   const Vector2D& p);

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

Conserved rotate_solve_rotate_back
(const RiemannSolver& rs,
 const Primitive& right,
 const Primitive& left,
 const double velocity,
 const Vector2D& n,
 const Vector2D& p);

#endif // SIMPLE_FLUX_CALCULATOR_HPP
