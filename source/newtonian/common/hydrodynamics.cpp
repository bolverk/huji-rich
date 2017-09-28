#include <cmath>
#include "hydrodynamics.hpp"

Primitive CalcPrimitive(double density, double pressure,
			Vector2D const& velocity,
			EquationOfState const& eos)
{
  return Primitive(density,
		   pressure,
		   velocity,
		   eos.dp2e(density, pressure),
		   eos.dp2c(density, pressure));
}

Primitive Conserved2Primitive
(Conserved const& c, EquationOfState const& eos)
{
  const double density = c.Mass;
  const Vector2D velocity = c.Momentum / c.Mass;
  const double energy = c.Energy / c.Mass -
	  0.5*ScalarProd(velocity, velocity);
  assert(energy>=0 &&
	 "Thermal energy is negative");
  /*
  if(energy<0)
  {
	  UniversalError eo("Negative thermal energy");
	  eo.AddEntry("x velocity",velocity.x);
	  eo.AddEntry("y velocity",velocity.y);
	  eo.AddEntry("density",density);
	  eo.AddEntry("Energy per volume",c.Energy);
	  throw eo;
  }
  */
  const double pressure = eos.de2p(density, energy);
  const double sound_speed = eos.dp2c(density, pressure);
  return Primitive(density,
		   pressure,
		   velocity,
		   energy,
		   sound_speed);
}

Primitive make_eos_consistent
(Primitive const& p,
 EquationOfState const& eos)
{
  return CalcPrimitive(p.Density,
		       p.Pressure,
		       p.Velocity,
		       eos);
}
