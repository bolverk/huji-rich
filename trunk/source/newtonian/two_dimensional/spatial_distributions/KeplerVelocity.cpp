#include "KeplerVelocity.hpp"

KeplerVelocity::KeplerVelocity(Direction dir,double Mass,Vector2D const&
	center):dir_(dir),M_(Mass),center_(center){}

double KeplerVelocity::operator()(Vector2D const& vec) const
{
	const double r=abs(center_-vec);
	if(dir_==xdir)
		return -sqrt(M_/pow(r,3))*(vec.y-center_.y);
	else
		return sqrt(M_/pow(r,3))*(vec.x-center_.x);
}
