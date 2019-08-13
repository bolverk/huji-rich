/*! \file Hllc3D.hpp
\brief HLLC riemann solver on an eulerian grid in 3D
\details This file is based on a code originally written by Omer Bromberg
\author Elad Steinberg
*/

#ifndef HLLC3DE_HPP
#define HLLC3DE_HPP 1

#include "RiemannSolver3D.hpp"

//! \brief HLLC Riemann solver for an Eulerian grid
class Hllc3DEnergy : public RiemannSolver3D
{
private:
	const double gamma_;
public:
	Hllc3DEnergy(double gamma = -1);

	Conserved3D operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
		EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const;
};

#endif //HLLC3DE_HPP
