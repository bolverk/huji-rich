/*! \file LagrangianHLLC3D.hpp
\brief HLLC riemann solver that keeps track of transfered internal energy
\author Elad Steinberg
*/

#ifndef LAGRANGIANHLLC3D_HPP
#define LAGRANGIANHLLC3D_HPP 1

#include "../three_dimensional/RiemannSolver3D.hpp"

//! \brief LagrangianHLLC Riemann solver for an Eulerian grid
class LagrangianHLLC3D : public RiemannSolver3D
{
private:
	const bool massflux_;
public:
	/*! \brief Class constructor
	\param massflux Whether to apply correction for mass flux
	*/
	explicit LagrangianHLLC3D(bool massflux = true);

	Conserved3D operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
		EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const;

	//! \brief Energy
	mutable double ws;
};

#endif //LAGRANGIANHLLC3D_HPP