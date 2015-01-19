/*! \file FirstOrderHydroFlux.hpp
  \brief Class for flux calculator for first order hydro. Assumes rigid walls as default.
  \author Elad Steinberg
 */

#ifndef FIRST_ORDER_HYDRO_FLUX_HPP
#define FIRST_ORDER_HYDRO_FLUX_HPP 1


#include "flux_calculator.hpp"
#include "../common/riemann_solver.hpp"
#include "../../misc/utils.hpp"

//! \brief First order flux calculator
class FirstOrderHydroFlux : public FluxCalculator
{
private:
	RiemannSolver const& rs_;
public:
	/*!
	\brief Class constructor
	\param rs The Riemann solver
	*/
	FirstOrderHydroFlux(RiemannSolver const& rs);

	//! \brief Class destructor
	~FirstOrderHydroFlux(void);

	vector<Conserved3D> operator()(const Tessellation3D& tess,
		const vector<ComputationalCell>& cells,const EquationOfState& eos,
		const vector<Vector3D>& point_velocities) const;


};

#endif //FIRST_ORDER_HYDRO_FLUX_HPP
