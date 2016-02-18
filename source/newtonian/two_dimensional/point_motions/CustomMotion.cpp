#include "CustomMotion.hpp"

CustomMotionCriteria::~CustomMotionCriteria(void){}

CustomMotion::CustomMotion(PointMotion const& otherpm, CustomMotionCriteria const& criteria) : pm_(otherpm),
criteria_(criteria){}

vector<Vector2D> CustomMotion::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
	double time, TracerStickerNames const& tracerstickernames) const
{
	return pm_.operator()(tess, cells, time,tracerstickernames);
}

vector<Vector2D> CustomMotion::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	double dt, vector<Vector2D> const& velocities, TracerStickerNames const& tracerstickernames)const
{
	vector<Vector2D> res =  pm_.ApplyFix(tess, cells, time, dt, velocities,tracerstickernames);
	for (size_t i = 0; i < res.size(); ++i)
	{
		if (criteria_.SatisfyCriteria(i, tess, cells, time,velocities,dt,tracerstickernames))
			res[i] = criteria_.CustomVelocityResult(i, tess, cells, time,velocities,dt,tracerstickernames);
	}
	return res;
}
