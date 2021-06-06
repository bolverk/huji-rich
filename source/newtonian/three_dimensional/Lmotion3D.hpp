#ifndef LMOTION3D_HPP
#define LMOTION3D_HPP 1
#include "point_motion_3d.hpp"
#include "LinearGauss3D.hpp"

//! \brief Point motion according to Lloyd's algorithm
class LMotion3D : public PointMotion3D
{
private:
	LinearGauss3D const& interp_;
	EquationOfState const& eos_;
	const double round_speed_,max_v_correction_;
public:
	/*! \brief Class constructor
	\param max_v_correction The maximum fractional correction to add to the velocity because of fix from previous time steps
	\param interp Interpolation scheme
	\param eos Equation of state
	\param roundspeed Speed of correction speed
	*/
	LMotion3D(LinearGauss3D const& interp, EquationOfState const& eos,double roundspeed=0.1,
		double max_v_correction=0.05);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, vector<Vector3D> &res) const override;

	void ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
		double dt, vector<Vector3D> &velocities)const override;
};
#endif //LMOTION3D_HPP
