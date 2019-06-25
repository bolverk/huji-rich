/*! \file Hllc3D.hpp
\brief HLLC riemann solver on an eulerian grid in 3D
\details This file is based on a code originally written by Omer Bromberg
\author Elad Steinberg
*/

#ifndef HLLC3D_HPP
#define HLLC3D_HPP 1

#include "RiemannSolver3D.hpp"

//! \brief HLLC Riemann solver for an Eulerian grid
class Hllc3D : public RiemannSolver3D
{
private:
	const double gamma_;
public:
	Hllc3D(double gamma = -1);

	Conserved3D operator()(ComputationalCell3D const& left,	ComputationalCell3D const& right,double velocity,
		EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const;
};

#endif //HLLC3D_HPP
