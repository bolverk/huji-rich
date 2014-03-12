#include "RotationForce.hpp"


RotationForce::RotationForce(Vector2D const& origin,Vector2D const& direction,
	int tracer_index):
origin_(origin),direction_(direction/abs(direction)),tracer_index_(tracer_index)
{}


RotationForce::~RotationForce(){}


Conserved RotationForce::Calculate
	(Tessellation const& tess,
	vector<Primitive> const& cells,
	 int point,vector<Conserved> const& /*fluxes*/,
	vector<Vector2D> const& /*point_velocity*/,
	 HydroBoundaryConditions const& /*hbc*/,
	vector<vector<double> > const & tracer,vector<double> & /*dtracer*/,
	double /*time*/,double /*dt*/)
{
	const Vector2D point_vec(tess.GetCellCM(point));
	const double r=ScalarProd(point_vec-origin_,direction_);
	const double mass=tess.GetVolume(point)*cells[point].Density;
	const Vector2D momentum=mass*tracer[point][tracer_index_]*
		tracer[point][tracer_index_]*direction_/(r*r*r);
	const double energy=ScalarProd(momentum,cells[point].Velocity);
	return Conserved(0,momentum,energy);
}
