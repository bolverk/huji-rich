/*! \file CustomMotion.hpp
\brief Abstract class for a custom motion of the mesh generating points
\author Elad Steinberg
*/

#ifndef CUSTOM_MOTION_HPP
#define CUSTOM_MOTION_HPP 1

#include "../point_motion.hpp"

//! \brief Class for checking if the criteria for custom motion is applied
class CustomMotionCriteria
{
public:
	virtual bool SatisfyCriteria(size_t index, Tessellation const& tess, vector<ComputationalCell> const& cells,
		double time)const = 0;
	virtual Vector2D CustomVelocityResult(size_t index, Tessellation const& tess, vector<ComputationalCell> const& cells,
		double time)const = 0;
	virtual ~CustomMotionCriteria(void);
};

/*! \brief Abstract class for custom motion of mesh generating points
*/
class CustomMotion : public PointMotion
{
public:
	CustomMotion(PointMotion const& otherpm, CustomMotionCriteria const& criteria);

	vector<Vector2D> operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
		double time) const;

	void ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
		double dt, vector<Vector2D> & velocities)const;
private:
	PointMotion const& pm_;
	CustomMotionCriteria const& criteria_;
};

#endif
