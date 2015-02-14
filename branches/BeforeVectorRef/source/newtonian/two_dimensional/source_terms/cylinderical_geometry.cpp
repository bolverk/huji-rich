#include "cylinderical_geometry.hpp"

CylindericalGeometry::CylindericalGeometry
(Vector2D const& origin,
 Vector2D const& direction):
  origin_(origin), direction_(direction) {}

double distance_from_line(Vector2D const& point,
			  Vector2D const& origin,
			  Vector2D const& direction)
{
  const double hypotenuse = abs(point - origin);
  const double side = fabs(Projection(point-origin,direction));
  return sqrt(pow(hypotenuse,2)-pow(side,2));
}

Vector2D cross_z(Vector2D const& v)
{
  return Vector2D(v.y,-v.x);
}

Conserved CylindericalGeometry::Calculate
(Tessellation const& tess,
 const PhysicalGeometry& /*pg*/,
 vector<Primitive> const& cells,
 int point,
 vector<Conserved> const& /*fluxes*/,
 vector<Vector2D> const& /*point_velocity*/,
 HydroBoundaryConditions const& /*hbc*/,
 vector<vector<double> > const &tracer,vector<double> &dtracer,
 vector<double> const& /*lengthes*/,
 double /*t*/,
 double /*dt*/)
{
  const double d = cells[static_cast<size_t>(point)].Density;
  const double vz = Projection(cells[static_cast<size_t>(point)].Velocity,
			       direction_);
  const Vector2D r_hat = cross_z(direction_);
  const double vr = Projection(cells[static_cast<size_t>(point)].Velocity,
			       r_hat);
  const double p = cells[static_cast<size_t>(point)].Pressure;
  const double e = cells[static_cast<size_t>(point)].Energy;
  const double r = distance_from_line(tess.GetCellCM(point),
				      origin_,
				      direction_);
  const double r_mom = d*pow(vr,2)/r;
  const double z_mom = d*vr*vz/r;
  const double volume = tess.GetVolume(point);
  if(!dtracer.empty())
	  dtracer=-d*volume*(vr/r)*tracer[static_cast<size_t>(point)];
  return -volume*Conserved
    (d*vr/r,
     r_mom*r_hat/abs(r_hat)+
     z_mom*direction_/abs(direction_),
     vr*(p+d*(e+0.5*pow(vr,2)+0.5*pow(vz,2)))/r);
}
