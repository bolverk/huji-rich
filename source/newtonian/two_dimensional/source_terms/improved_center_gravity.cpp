#include "improved_center_gravity.hpp"

ImprovedCenterGravity::ImprovedCenterGravity(double M,double Rmin,Vector2D center):
M_(M),Rmin_(Rmin),_center(center){}

Vector2D ImprovedCenterGravity::Calculate
(Tessellation const& tess,
 vector<Primitive> const& /*cells*/,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const& /*tracers*/,
 double /*time*/,
 double /*dt*/)
{
	Vector2D pos(tess.GetCellCM(point)-_center);
	double r=abs(pos);
	return Vector2D((-1)*pos*M_/(r*r*r+Rmin_*Rmin_*Rmin_));
}

void ImprovedCenterGravity::set_center(Vector2D const& new_center)
{
  _center = new_center;
}

Vector2D const& ImprovedCenterGravity::get_center(void) const
{
  return _center;
}