#ifndef CENTRIPETAL_HPP
#define CENTRIPETAL_HPP 1

#include "../source_terms/ConservativeForce.hpp"

class Centripetal: public Acceleration
{
public:

  Centripetal(double angular_velocity,
	      Vector2D const& center);

  Vector2D Calculate(Tessellation const* tess,
		     vector<Primitive> const& /*cells*/,
		     int point,
		     vector<Conserved> const& /*fluxes*/,
		     vector<Vector2D> const& /*point_velocity*/,
		     HydroBoundaryConditions const* /*hbc*/,
		     double /*time*/,
		     double /*dt*/);

private:
  const double omega_;
  const Vector2D center_;
};

#endif // CENTRIPETAL_HPP
