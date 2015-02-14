#include "RotationForce.hpp"

RotationForce::RotationForce(Vector2D const& origin,Vector2D const& direction,
	int tracer_index):
origin_(origin),direction_(direction/abs(direction)),tracer_index_(tracer_index)
{}

RotationForce::~RotationForce(){}

Conserved RotationForce::Calculate
	(Tessellation const& tess,
	 const PhysicalGeometry& /*pg*/,
	vector<Primitive> const& cells,
	 int point,vector<Conserved> const& /*fluxes*/,
	vector<Vector2D> const& /*point_velocity*/,
	 HydroBoundaryConditions const& /*hbc*/,
	vector<vector<double> > const & tracer,vector<double> & /*dtracer*/,
	vector<double> const& /*lengthes*/,
	double /*time*/,double /*dt*/)
{
	const Vector2D& point_vec = tess.GetCellCM(point);
	const double r=ScalarProd(point_vec-origin_,direction_);
	const double mass=tess.GetVolume(point)*cells[static_cast<size_t>(point)].Density;
	const Vector2D momentum=mass*tracer[static_cast<size_t>(point)][static_cast<size_t>(tracer_index_)]*
	  tracer[static_cast<size_t>(point)][static_cast<size_t>(tracer_index_)]*direction_/(r*r*r);
	const double energy=ScalarProd(momentum,cells[static_cast<size_t>(point)].Velocity);
	return Conserved(0,momentum,energy);
}
