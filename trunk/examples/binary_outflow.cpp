#include "binary_outflow.hpp"

BinaryOutflowVelocity::BinaryOutflowVelocity(double R,double val,Vector2D loc,Dir axis):
R_(R),val_(val),loc_(loc),axis_(axis)
{}

double BinaryOutflowVelocity::operator()(Vector2D const& point) const
{
	const double d=point.distance(loc_);
	if(d<R_&&d>0)
	{
		Vector2D direction=point-loc_;
		direction=direction/abs(direction);
		return val_*((axis_==Xaxis) ? direction.x : direction.y);
	}
	else
		return 0;
}
