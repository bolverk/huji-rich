#include <cmath>
#include "hydrodynamic_variables.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"

Primitive::Primitive(void):
  Density(0),
  Pressure(0),
  Velocity(0,0),
  Energy(0),
  SoundSpeed(0) {}

Primitive::Primitive(const Primitive& other):
  Density(other.Density),
  Pressure(other.Pressure),
  Velocity(other.Velocity),
  Energy(other.Energy),
  SoundSpeed(other.SoundSpeed) {}

Primitive::Primitive(double density_i,
		     double pressure_i,
		     Vector2D const& velocity_i,
		     double energy_i,
		     double sound_speed_i):
  Density(density_i),
  Pressure(pressure_i),
  Velocity(velocity_i),
  Energy(energy_i),
  SoundSpeed(sound_speed_i) {}

int Primitive::GetVarNo(void) const
{
  return 6;
}

Primitive& Primitive::operator=(Primitive const&other)
{
  if (this == &other)
    return *this;
  Density=other.Density;
  Energy=other.Energy;
  Pressure=other.Pressure;
  SoundSpeed=other.SoundSpeed;
  Velocity=other.Velocity;
  return *this;
}

bool primitive_has_nan(Primitive const& p)
{
  return is_nan(p.Density)||
    is_nan(p.Pressure)||
    is_nan(p.Velocity.x)||
    is_nan(p.Velocity.y)||
    is_nan(p.Energy)||
    is_nan(p.SoundSpeed);
}

Conserved::Conserved(double mass,
		     Vector2D const& momentum,
		     double energy):
  Mass(mass), Momentum(momentum),
  Energy(energy) {}

Conserved::Conserved(const Conserved& other):
  Mass(other.Mass),
  Momentum(other.Momentum),
  Energy(other.Energy) {}

Conserved& Conserved::operator=(Conserved const&other)
{
  if (this == &other)
    return *this;
  Energy=other.Energy;
  Mass=other.Mass;
  Momentum=other.Momentum;
  return *this;
}

Primitive& Primitive::operator+=(Primitive const& p)
{
  Density += p.Density;
  Pressure += p.Pressure;
  Energy += p.Energy;
  SoundSpeed += SoundSpeed;
  Velocity.x += p.Velocity.x;
  Velocity.y += p.Velocity.y;
  return *this;
}

double& Primitive::operator[](int index)
{
  if(index==0)
    return Density;
  else if(index==1)
    return Pressure;
  else if(index==2)
    return Energy;
  else if(index==3)
    return SoundSpeed;
  else if(index==4)
    return Velocity.x;
  else if(index==5)
    return Velocity.y;
  else
    throw UniversalError("Invalid index in Primitive::operator[]");
}

double Primitive::operator[](int index) const
{
  if(index==0)
    return Density;
  else if(index==1)
    return Pressure;
  else if(index==2)
    return Energy;
  else if(index==3)
    return SoundSpeed;
  else if(index==4)
    return Velocity.x;
  else if(index==5)
    return Velocity.y;
  else
    throw UniversalError("Invalid index in Primitive::operator[]");
}

namespace{
  Primitive binary_op_primitive(Primitive const& p1,
				Primitive const& p2,
				BinaryOperation<double> const& op)
  {
    return Primitive(op(p1.Density,p2.Density),
		     op(p1.Pressure,p2.Pressure),
		     Vector2D(op(p1.Velocity.x,p2.Velocity.x),
			      op(p1.Velocity.y,p2.Velocity.y)),
		     op(p1.Energy,p2.Energy),
		     op(p1.SoundSpeed,p2.SoundSpeed));
  }

  Primitive primitive_scalar_op(Primitive const& p,
				double c,
				BinaryOperation<double> const& op)
  {
    return Primitive(op(p.Density,c),
		     op(p.Pressure,c),
		     Vector2D(op(p.Velocity.x,c),
			      op(p.Velocity.y,c)),
		     op(p.Energy,c),
		     op(p.SoundSpeed,c));
  }
}

Primitive operator+(Primitive const& p1,
		    Primitive const& p2)
{
  class plus: public BinaryOperation<double>
  {
    double operator()(double const& x, double const& y) const
    {
      return x+y;
    }
  } op;
  return binary_op_primitive(p1, p2, op);
}

Primitive operator-(Primitive const& p1,
		    Primitive const& p2)
{
  class minus: public BinaryOperation<double>
  {
  public:
    double operator()(double const& x, double const& y) const
    {
      return x-y;
    }
  } op;

  return binary_op_primitive(p1,p2,op);
}

Primitive operator/(Primitive const& p,
		    double c)
{
  class divide: public BinaryOperation<double>
  {
  public:
    double operator()(double const& x,double const& y) const
    {
      return x/y;
    }
  } op;
  return primitive_scalar_op(p,c,op);
}

Primitive operator*(Primitive const& p,
		    double s)
{
  class multiply: public BinaryOperation<double>
  {
  public:
    double operator()(double const& x, double const& y) const
    {
      return x*y;
    }
  } op;
  return primitive_scalar_op(p,s,op);
}

Primitive operator*(double s,
		    Primitive const& p)
{
  return p*s;
}

Conserved::Conserved(void):
  Mass(0), Momentum(Vector2D(0,0)), Energy(0) {}

Conserved operator+(Conserved const& c1, Conserved const& c2)
{
  return Conserved(c1.Mass+c2.Mass,
		   c1.Momentum+c2.Momentum,
		   c1.Energy+c2.Energy);
}

Conserved operator-(Conserved const& c1, Conserved const& c2)
{
  return Conserved(c1.Mass-c2.Mass,
		   c1.Momentum-c2.Momentum,
		   c1.Energy-c2.Energy);
}

Conserved operator*(double d, Conserved const& c)
{
  return Conserved(d*c.Mass,
		   d*c.Momentum,
		   d*c.Energy);
}

Conserved operator/(Conserved const& c, double d)
{
  return Conserved(c.Mass/d,
		   c.Momentum/d,
		   c.Energy/d);
}

Conserved& Conserved::operator+=(Conserved const& c)
{
  Mass += c.Mass;
  Momentum += c.Momentum;
  Energy += c.Energy;
  return *this;
}

Conserved& Conserved::operator-=(Conserved const& c)
{
  Mass -= c.Mass;
  Momentum -= c.Momentum;
  Energy -= c.Energy;
  return *this;
}

double TotalEnergyDensity(Primitive const& p)
{
  return p.Density*(0.5*pow(abs(p.Velocity),2)+
		    p.Energy);
}

Conserved Primitive2Conserved(Primitive const& p)
{
  return Conserved(p.Density,
		   p.Density*p.Velocity,
		   TotalEnergyDensity(p));
}

Conserved Primitive2Conserved(Primitive const& p,double vol)
{
  return Conserved(vol*p.Density,
		   vol*p.Density*p.Velocity,
		   vol*TotalEnergyDensity(p));
}

Conserved Primitive2Flux(Primitive const& p,
			 Vector2D const& n)
{
  const Vector2D nn = n/abs(n);
  return Conserved(p.Density*ScalarProd(p.Velocity, nn),
		   p.Pressure*nn +
		   p.Density*ScalarProd(p.Velocity, nn)*p.Velocity,
		   (TotalEnergyDensity(p)+p.Pressure)*
		   ScalarProd(p.Velocity, nn));
}
