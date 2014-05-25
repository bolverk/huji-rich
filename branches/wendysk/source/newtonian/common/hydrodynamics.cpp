#include <cmath>
#include "hydrodynamics.hpp"

Primitive CalcPrimitive(double density, double pressure,
			Vector2D const& velocity, 
			EquationOfState const& eos)
{
  Primitive res;
  res.Density = density;
  res.Pressure = pressure;
  res.Velocity = velocity;
  res.Energy = eos.dp2e(res.Density, res.Pressure);
  res.SoundSpeed = eos.dp2c(res.Density, res.Pressure);
  return res;
}

Primitive Conserved2Primitive
(Conserved const& c, EquationOfState const& eos)
{
  Primitive res;
  res.Density = c.Mass;
  res.Velocity = c.Momentum / c.Mass;
  res.Energy = c.Energy/c.Mass - 0.5*pow(abs(res.Velocity),2);
  res.Pressure = eos.de2p(res.Density, res.Energy);
  res.SoundSpeed = eos.dp2c(res.Density, res.Pressure);
  return res;
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
