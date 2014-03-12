#include "CenterGravity.hpp"

CenterGravity::CenterGravity(double M,double Rmin,Vector2D center):
M_(M),Rmin_(Rmin),_center(center){}

Vector2D CenterGravity::Calculate
(Tessellation const& tess,
 vector<Primitive> const& /*cells*/,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 double /*time*/,
 double /*dt*/)
{
	Vector2D pos(tess.GetCellCM(point)-_center);
	double r=abs(pos);
	return Vector2D((-1)*pos*M_/(r*r*r+Rmin_*Rmin_*Rmin_));
}
