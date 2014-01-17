#include <cmath>
#include "hydrodynamic_variables.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"

Primitive::Primitive(void):
  Density(0), Pressure(0),
   Energy(0),
  Velocity(0,0), SoundSpeed(0) {}

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

//! \brief Operation on two numbers that reeturns a number
class BinaryOperation
{
public:

  /*! \brief Abstract class for binary operation
    \param x Left argument
    \param y Right argument
    \return Result of the binary operation
   */
  virtual double Eval(double x, double y) const = 0;

  virtual ~BinaryOperation(void) {}
};

namespace{
  Primitive binary_op_primitive(Primitive const& p1,
				Primitive const& p2,
				BinaryOperation const& op)
  {
    Primitive res;
    for(int i = 0;i<res.GetVarNo();i++)
      res[i] = op.Eval(p1[i],p2[i]);
    return res;
  }

  Primitive primitive_scalar_op(Primitive const& p,
				double c,
				BinaryOperation const& op)
  {
    Primitive res;
    for(int i=0;i<res.GetVarNo();i++)
      res[i] = op.Eval(p[i],c);
    return res;
  }
}

Primitive operator+(Primitive const& p1,
		   Primitive const& p2)
{
  class plus: public BinaryOperation
  {
    double Eval(double x, double y) const
    {
      return x+y;
    }
  } op;
  return binary_op_primitive(p1, p2, op);
}		  

Primitive operator-(Primitive const& p1,
		    Primitive const& p2)
{
  class minus: public BinaryOperation
  {
  public:
    double Eval(double x, double y) const
    {
      return x-y;
    }
  } op;

  return binary_op_primitive(p1,p2,op);
}

Primitive operator/(Primitive const& p,
		    double c)
{
  class divide: public BinaryOperation
  {
  public:
    double Eval(double x,double y) const
    {
      return x/y;
    }
  } op;
  return primitive_scalar_op(p,c,op);
}

Primitive operator*(Primitive const& p,
		    double s)
{
  class multiply: public BinaryOperation
  {
  public:
    double Eval(double x, double y) const
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
  Conserved res;
  res.Mass = c1.Mass + c2.Mass;
  res.Momentum = c1.Momentum + c2.Momentum;
  res.Energy = c1.Energy + c2.Energy;
  return res;
}

Conserved operator-(Conserved const& c1, Conserved const& c2)
{
  Conserved res;
  res.Mass = c1.Mass - c2.Mass;
  res.Momentum = c1.Momentum - c2.Momentum;
  res.Energy = c1.Energy - c2.Energy;
  return res;
}

Conserved operator*(double d, Conserved const& c)
{
  Conserved res;
  res.Mass = d * c.Mass;
  res.Momentum = d * c.Momentum;
  res.Energy = d * c.Energy;
  return res;
}

Conserved operator/(Conserved const& c, double d)
{
  Conserved res;
  res.Mass = c.Mass / d;
  res.Momentum = c.Momentum / d;
  res.Energy = c.Energy / d;
  return res;
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
  Conserved res;
  res.Mass = p.Density;
  res.Momentum = p.Density*p.Velocity;
  res.Energy = TotalEnergyDensity(p);
  return res;
}

Conserved Primitive2Conserved(Primitive const& p,double vol)
{
  Conserved res;
  res.Mass = vol*p.Density;
  res.Momentum = vol*p.Density*p.Velocity;
  res.Energy = vol*TotalEnergyDensity(p);
  return res;
}


Conserved Primitive2Flux(Primitive const& p, 
			 Vector2D const& n)
{
  Vector2D nn = n/abs(n);
  Conserved res;
  res.Mass = p.Density*ScalarProd(p.Velocity, nn);
  res.Momentum = p.Pressure*nn + 
    p.Density*ScalarProd(p.Velocity, nn)*p.Velocity;
  res.Energy = (TotalEnergyDensity(p)+p.Pressure)*
    ScalarProd(p.Velocity, nn);
  return res;
}
