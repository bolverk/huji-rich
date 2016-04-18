/*! \file simple_flux_calculator.hpp
  \author Almog Yalinewich
  \brief Simple flux calculator
 */

#ifndef SIMPLE_FLUX_CALCULATOR_HPP
#define SIMPLE_FLUX_CALCULATOR_HPP 1

#include "flux_calculator_2d.hpp"
#include "../common/riemann_solver.hpp"

 /*! \brief Converts computational cell to primitive variables
   \param cell Computational cell
   \param eos Equation of state
   \param tracerstickernames The names of the tracers and stickers
   \return Primitive variable
  */
Primitive convert_to_primitive(const ComputationalCell& cell,
	const EquationOfState& eos, TracerStickerNames const& tracerstickernames);

/*! \brief Reflects velocity about axis
  \param p Primitive variables
  \param axis Reflection axis
  \return Primitive variables with reflected velocity
 */
Primitive reflect(const Primitive& p,
	const Vector2D& axis);

/*! \brief Remove parallel component of a vector
  \param v Vector
  \param p Parallel direction
  \return v with parallel component removed
 */
Vector2D remove_parallel_component(const Vector2D& v,
	const Vector2D& p);

//! \brief Simple flux calculator
class SimpleFluxCalculator : public FluxCalculator
{
public:

	/*! \brief Class constructor
	  \param rs Riemann solver
	 */
	explicit SimpleFluxCalculator(const RiemannSolver& rs);

	vector<Extensive> operator()
		(const Tessellation& tess,
			const vector<Vector2D>& edge_velocities,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& extensives,
			const CacheData& cd,
			const EquationOfState& eos,
			const double time,
			const double dt,
			TracerStickerNames const& tracerstickernames) const;

private:
	const RiemannSolver& rs_;

	Conserved calcHydroFlux(const Tessellation& tess,
		const vector<Vector2D>& edge_velocities,
		const vector<ComputationalCell>& cells,
		const EquationOfState& eos,
		const size_t i,
		TracerStickerNames const& tracerstickernames) const;
};

/*! \brief Rotates, solve riemann problem and rotates results back
  \param rs Riemann solver
  \param left Primitive variables on the left side
  \param right Primitive variables on the right side
  \param velocity Velocity of the interface
  \param n Normal direction
  \param p Parallel direction
  \return Hydrodynamic fluxes
 */
Conserved rotate_solve_rotate_back
(const RiemannSolver& rs,
	const Primitive& left,
	const Primitive& right,
	const double velocity,
	const Vector2D& n,
	const Vector2D& p);

#endif // SIMPLE_FLUX_CALCULATOR_HPP
