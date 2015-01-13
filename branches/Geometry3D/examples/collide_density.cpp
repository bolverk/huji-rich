#include "collide_density.hpp"

CollideDensity::CollideDensity(double density1,double density2,double background,
	double x1,double y1,double x2,double y2,double r1,double r2):
density1_(density1),density2_(density2),background_(background),x1_(x1),
	y1_(y1),x2_(x2),y2_(y2),r1_(r1),r2_(r2){}

double CollideDensity::operator()(Vector2D const& point) const
{
	if(point.distance(Vector2D(x1_,y1_))<r1_)
		return density1_;
	if(point.distance(Vector2D(x2_,y2_))<r2_)
		return density2_;
	return background_;
}
