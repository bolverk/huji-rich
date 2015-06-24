#include "CenterGravity.hpp"

CenterGravity::CenterGravity(double M, double Rmin, double SoftStart,Vector2D center,double omega) :
M_(M),Rmin_(Rmin),softlength_(SoftStart),_center(center),omega_(omega){}

Vector2D CenterGravity::Calculate
(Tessellation const& tess,
 vector<Primitive> const& /*cells*/,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const& /*tracers*/,
 double time,
 double /*dt*/)
{
	Vector2D pos(tess.GetCellCM(point) - (abs(_center)*Vector2D(cos(omega_*time), sin(omega_*time))));
	//Vector2D pos(tess.GetCellCM(point)-_center);
	double r=abs(pos);
	if (r<softlength_)
		return Vector2D((-1)*pos*M_*pow(r*r+Rmin_*Rmin_,-1.5));
	else
		return Vector2D((-1)*pos*M_ / (r*r*r));
}
