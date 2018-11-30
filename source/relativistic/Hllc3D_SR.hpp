/*! \file Hllc3D_SR.hpp
\brief HLLC riemann solver in 3D with SR
\author Elad Steinberg
*/

#ifndef HLLC3D_SR_HPP
#define HLLC3D_SR_HPP 1

#include "../newtonian/three_dimensional/RiemannSolver3D.hpp"

//! \brief HLLC Riemann solver for an Eulerian grid
class Hllc3D_SR : public RiemannSolver3D
{
private:
	mutable ComputationalCell3D local_left_, local_right_;
public:
	Hllc3D_SR();

	~Hllc3D_SR();

	Conserved3D operator()(ComputationalCell3D const& left, ComputationalCell3D const& right, double velocity,
		EquationOfState const& eos, TracerStickerNames const& tsn, Vector3D const& normaldir) const;
};

#endif //HLLC3D_HPP
