/*! \file RoundCElls.hpp
\author Elad Steinberg
\brief Correction to point velocities that keeps cells round
*/

#ifndef ROUND_CELLS3D_HPP
#define ROUND_CELLS3D_HPP 1

#include "point_motion_3d.hpp"
#include "../common/equation_of_state.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"

//! \brief Correction to point velocities that keeps cells round
//! \details Based on equation 63 in the Arepo paper
class RoundCells3D : public PointMotion3D
{
public:

	/*! \brief Class constructor
	\param pm Base point motion
	\param eos Equation of state
	\param outer The outer boudnary conditions, used for preventing points from getting outside the box. If periodic then no fix is applied
	\param chi chi parameter in equation 63
	\param eta eta parameter in equation 63
	\param cold Switch for cold flows
	*/
	RoundCells3D(const PointMotion3D& pm, const EquationOfState& eos,Vector3D const& ll,Vector3D const& ur,
		double chi = 0.25, double eta = 0.02, bool cold = false);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const;

	void ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
		double dt, vector<Vector3D> &velocities, TracerStickerNames const& tracerstickernames)const;
private:

	void calc_dw(Vector3D &velocty, size_t i, const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		TracerStickerNames const& tracerstickernames,vector<Vector3D> & velocities) const;

	void calc_dw(Vector3D &velocty, size_t i, const Tessellation3D& tess, double dt, vector<ComputationalCell3D> const& cells,
		TracerStickerNames const& tracerstickernames, vector<Vector3D> & velocities)const;

	const PointMotion3D& pm_;
	const EquationOfState& eos_;
	const Vector3D ll_, ur_;
	const double chi_;
	const double eta_;
	const bool cold_;
};

#endif // ROUND_CELLS3D_HPP
