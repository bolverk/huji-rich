#include "SimpleGravity.hpp"

SimpleGravity::SimpleGravity(double M, double Rmin, double SoftStart, Vector2D center,double omega) :
M_(M), Rmin_(Rmin), softlength_(SoftStart), _center(center),omega_(omega),first_time_(true),dt_(-1){}

Conserved SimpleGravity::Calculate
(Tessellation const& tess,
const PhysicalGeometry& pg,
vector<Primitive> const& cells,
int point,
vector<Conserved> const& /*fluxes*/,
vector<Vector2D> const& /*point_velocity*/,
HydroBoundaryConditions const& /*hbc*/,
vector<vector<double> > const& /*tracers*/,
vector<double>& /*dtracer*/, vector<double> const& /*lengthes*/,
double time,
double dt)
{
	Vector2D pos(tess.GetCellCM(point) - (abs(_center)*Vector2D(cos(omega_*time),sin(omega_*time))));
	double r = abs(pos);
	double m = tess.GetVolume(point)*cells[point].Density;
	Vector2D acc;
	if (r<softlength_)
		acc=(-1)*pos*M_ / (r*r*r + Rmin_*Rmin_*r);
	else
		acc=((-1)*pos*M_ / (r*r*r));
	if (first_time_)
	{
		dt_ = sqrt(tess.GetWidth(point) / abs(acc));
		first_time_ = false;
		time_ = time;
	}
	else
	{
		if (time_ < time)
		{
			dt_ = sqrt(tess.GetWidth(point) / abs(acc));
			time_ = time;
		}
		else
			dt_ = std::min(dt_, sqrt(tess.GetWidth(point) / abs(acc)));
	}
	return Conserved(0, m*acc, m*ScalarProd(acc, cells[point].Velocity));
}

double SimpleGravity::GetTimeStep(void) const
{
	return dt_;
}