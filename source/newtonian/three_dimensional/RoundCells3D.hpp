/*! \file RoundCells3D.hpp
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
	\param chi chi parameter in equation 63
	\param eta eta parameter in equation 63
	\param cold Switch for cold flows
	\param min_dw The minimum magnitude for the cell' velocity scale
	\param no_move List of stickers not to move the cells for
	\param dt_speed The speed in units of radius/dt
	*/
	RoundCells3D(const PointMotion3D& pm, const EquationOfState& eos,
		double chi = 0.25, double eta = 0.02, bool cold = false,double min_dw=0,double dt_speed=0.01,
		const vector<std::string>& no_move=vector<std::string>());

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, vector<Vector3D> &res) const override;

	void ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
		double dt, vector<Vector3D> &velocities)const override;
private:

	void calc_dw(Vector3D &velocty, size_t i, const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Vector3D> & velocities, vector<char> const& nomove) const;

	void calc_dw(Vector3D &velocty, size_t i, const Tessellation3D& tess, double dt, vector<ComputationalCell3D> const& cells,
		vector<Vector3D> & velocities, vector<char> const& nomove)const;

	const PointMotion3D& pm_;
	const EquationOfState& eos_;
	const double chi_;
	const double eta_;
	const bool cold_;
	const double min_dw_;
	const double dt_speed_;
	const vector<std::string> no_move_;
};

#endif // ROUND_CELLS3D_HPP
