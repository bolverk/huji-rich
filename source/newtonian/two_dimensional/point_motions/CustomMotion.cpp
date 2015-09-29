#include "CustomMotion.hpp"

CustomMotionCriteria::~CustomMotionCriteria(void){}

CustomMotion::CustomMotion(PointMotion const& otherpm, CustomMotionCriteria const& criteria) : pm_(otherpm),
criteria_(criteria){}

vector<Vector2D> CustomMotion::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
	double time) const
{
	return pm_.operator()(tess, cells, time);
}

void CustomMotion::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	double dt, vector<Vector2D> & velocities)const
{
	pm_.ApplyFix(tess, cells, time, dt, velocities);
	for (size_t i = 0; i < velocities.size(); ++i)
	{
		if (criteria_.SatisfyCriteria(i, tess, cells, time))
			velocities[i] = criteria_.CustomVelocityResult(i, tess, cells, time);
	}
}