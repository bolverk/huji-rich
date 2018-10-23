#ifndef LMOTION3D_HPP
#define LMOTION3D_HPP 1
#include "point_motion_3d.hpp"
#include "LinearGauss3D.hpp"

class LMotion3D : public PointMotion3D
{
private:
	LinearGauss3D const& interp_;
	EquationOfState const& eos_;
	const double round_speed_,max_v_correction_;
public:
	/*!
	\param max_v_correction The maximum fractional correction to add to the velocity because of fix from previous time steps
	*/
	LMotion3D(LinearGauss3D const& interp, EquationOfState const& eos,double roundspeed=0.1,
		double max_v_correction=0.05);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const;

	void ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
		double dt, vector<Vector3D> &velocities, TracerStickerNames const& tracerstickernames)const;
};
#endif //LMOTION3D_HPP