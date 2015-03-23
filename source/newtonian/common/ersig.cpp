#include <algorithm>
#include <cmath>
#include "ersig.hpp"
#include "ideal_gas.hpp"
#include "hydrodynamics.hpp"
#include "../../misc/bisection.hpp"
#include "../../misc/universal_error.hpp"

using namespace std;

namespace {
  double phugo(double v2,
	       double d1,
	       double p1,
	       double g)
  {
    return (4*p1 + d1*(1 + g)*pow(v2,2) + sqrt(d1)*v2*
	    sqrt(16*g*p1 + d1*pow(1 + g,2)*pow(v2,2)))/4.;
  }

  double vhugo(double p2,
	       double d1,
	       double p1,
	       double g)
  {
    return sqrt(2.)*(p2-p1)/sqrt(d1*(p2*(g+1)+p1*(g-1)));
  }

  double visen(double p2,
	       double d1,
	       double p1,
	       double g)
  {
    return (2*sqrt(g)*pow(p1,1/g/2)*
	    (pow(p2,(g-1)/g/2)-pow(p1,(g-1)/g/2)))/
      (sqrt(d1)*(g-1));
  }

  double vhydro(double p2,
		double d1,
		double p1,
		double g)
  {
    if(p2>p1)
      return vhugo(p2,d1,p1,g);
    else
      return visen(p2,d1,p1,g);
  }

  double left_velocity(double p,
		       Primitive const& upstream,
		       double g)
  {
    const double d1 = upstream.Density;
    const double p1 = upstream.Pressure;
    const double v1 = upstream.Velocity.x;
    return v1 - vhydro(p,d1,p1,g);
  }

  double right_velocity(double p,
			Primitive const& upstream,
			double g)
  {
    const double d1 = upstream.Density;
    const double p1 = upstream.Pressure;
    const double v1 = upstream.Velocity.x;
    return v1 + vhydro(p,d1,p1,g);
  }

  double eval_trans_eqn(double p,
			Primitive const& left,
			Primitive const& right,
			double g)
  {
    return right_velocity(p,right,g) - left_velocity(p,left,g);
  }

  //! \brief Transcendental equation for the pressure at the contact discontinuity
  class TransEqn: public Func1Var
  {
  public:

    /*! \brief Class constructor
      \param left Hydrodynamic variables on the left side
      \param right Hydrodynamic variables on the right side
      \param g Adiabatic index
    */
    TransEqn(Primitive const& left,
	     Primitive const& right,
	     double g):
      left_(left),
      right_(right),
      g_(g) {}

    /*! \brief Evaluates the equation
      \param p Pressure
      \return Value of the equation
    */
    double eval(double p) const
    {
      return eval_trans_eqn(p,left_,right_,g_);
    }

  private:

    const Primitive left_;
    const Primitive right_;
    const double g_;
  };

  double ps_two_rarefactions(Primitive const& left,
			     Primitive const& right,
			     double g)
  {
    const double dl = left.Density;
    const double pl = left.Pressure;
    const double vl = left.Velocity.x;
    const double dr = right.Density;
    const double pr = right.Pressure;
    const double vr = right.Velocity.x;
    const double temp =
      (2*(sqrt(g*dr*pl)+sqrt(g*dl*pr))+
       sqrt(dl*dr)*(g-1)*(vl-vr))/
      (2*(sqrt(g*dr)*pow(pl,1/g/2)+
	  sqrt(g*dl)*pow(pr,1/g/2)));
    if(temp<0){
      throw UniversalError("Vacuum formed");
    }
    const double res = pow(temp,2*g/(g-1));
    return res;
  }
}

namespace{
  class PVContact
  {
  public:

    /*! \brief Class constructor
      \param p Pressure at the contact discontinuity
      \param v Velocity at the contact discontinuity
    */
    PVContact(double p, double v):
      pressure(p), velocity(v) {}

    //! \brief Pressure at the contact discontinuity
    const double pressure;

    //! \brief Velocity at the contact discontinuity
    const double velocity;
  };
}

namespace{
  PVContact calc_contact_pv(Primitive const& left,
			    Primitive const& right,
			    double g)
  {
    using namespace std;

    const double tol = 1e-6;

    const double vlpr = left_velocity(right.Pressure, left, g);
    const double vrpl = right_velocity(left.Pressure, right, g);
    double ps = 0;
    if(min(left.Velocity.x,vlpr)>
       max(right.Velocity.x,vrpl)){
      // Two shocks
      const double plvr = phugo(left.Velocity.x-right.Velocity.x,
				left.Density,left.Pressure,g);
      const double prvl = phugo(left.Velocity.x-right.Velocity.x,
				right.Density,right.Pressure,g);
      const TransEqn eqn(left,right,g);
      try{
	ps = bisection(eqn,max(left.Pressure,right.Pressure),
		       2*max(plvr,prvl),tol);
      }
      catch(UniversalError& eo){
	eo.AddEntry("caught in two shock section of ersig",0);
	eo.AddEntry("left density",left.Density);
	eo.AddEntry("right density",right.Density);
	eo.AddEntry("left pressure",left.Pressure);
	eo.AddEntry("right pressure",right.Pressure);
	eo.AddEntry("left x velocity",left.Velocity.x);
	eo.AddEntry("right x velocity",right.Velocity.y);
	eo.AddEntry("rethrown",0);
	throw;
      }
    }
    else if(max(left.Velocity.x,vlpr)<
	    min(right.Velocity.x,vrpl)){
      // Two rarefactions
      ps = ps_two_rarefactions(left, right, g);
    }
    else{
      // Shock rarefaction
      const TransEqn eqn(left, right, g);
      ps = bisection(eqn,
		     min(left.Pressure, right.Pressure),
		     max(left.Pressure, right.Pressure),
		     tol);
    }
    const double vs = 0.5*(left_velocity(ps,left,g)+right_velocity(ps,right,g));
    return PVContact(ps,vs);
  }
}

// Self similar spatial profiles

// Shocks

namespace {
  double shock_speed(double p2,
		     double d1,
		     double p1,
		     double g)
  {
    return sqrt(g*(p2+p1)+p2-p1)/sqrt(2*d1);
  }

  double dhugo(double p2,
	       double d1,
	       double p1,
	       double g)
  {
    return (d1*((-1 + g)*p1 + (1 + g)*p2))/((1 + g)*p1 + (-1 + g)*p2);
  }

  Primitive move_on_hugoniot(double ps,
			     Primitive const& upstream,
			     double g)
  {
    const double ds  = dhugo(ps,
			     upstream.Density,
			     upstream.Pressure, g);
    const Vector2D vs(vhugo(ps,
			    upstream.Density,
			    upstream.Pressure,g),
		      upstream.Velocity.y);
    return CalcPrimitive(ds,ps, vs, IdealGas(g));
  }

  class HydroProf
  {
  public:

    virtual Primitive getHydroVars(double v) const = 0;

    virtual ~HydroProf(void) {}
  };

  Primitive x_rest_frame(Primitive const& source)
  {
    Primitive res =source;
    res.Velocity.x = 0;
    return res;
  }
}

//! \brief Spatial profile in case of shock wave
namespace {
  class ShockProf
  {
  public:

    ShockProf
    (double p2,
     Primitive const& upstream,
     double g):
      vs_(shock_speed(p2,
		      upstream.Density,
		      upstream.Pressure,g)),
      upstream_(x_rest_frame(upstream)),
      downstream_
      (move_on_hugoniot(p2,upstream,g)) {}

    Primitive getHydroVars(double v) const
    {
      if(v>vs_)
	return upstream_;
      else
	return downstream_;
    }

  private:

    const double vs_;
    const Primitive upstream_;
    const Primitive downstream_;
  };
}

/*class ShockProf: public HydroProf
  {
  public:

  ShockProf(double p2,
  Primitive const& upstream,
  double g):
  vs_(shock_speed(p2,
  upstream.Density,
  upstream.Pressure,g)),
  upstream_(x_rest_frame(upstream)),
  downstream_
  (move_on_hugoniot(p2,upstream,g)) {}

  Primitive getHydroVars(double v) const
  {
  if(v>vs_)
  return upstream_;
  else
  return downstream_;
  }

  private:

  const double vs_;
  const Primitive upstream_;
  const Primitive downstream_;
  };*/

// Isentropes

namespace {
  double disen(double p2,
	       double d1,
	       double p1,
	       double g)
  {
    return d1*pow(p2/p1,1/g);
  }

  Primitive move_on_isentrope(double pr,
			      Primitive const& upstream,
			      double g)
  {
    const double dr = disen(pr,
			    upstream.Density,
			    upstream.Pressure,g);
    const Vector2D vr(visen(pr,
			    upstream.Density,
			    upstream.Pressure,g),
		      upstream.Velocity.y);
    return CalcPrimitive(dr,pr,vr,IdealGas(g));
  }
}

namespace {
  class IsenProf
  {
  public:

    IsenProf
    (double p2,
     Primitive const& upstream,
     double g):
      g_(g),
      upstream_(x_rest_frame(upstream)),
      downstream_
      (move_on_isentrope(p2, upstream,g)) {}

    Primitive getHydroVars(double v) const
    {
      if(v>upstream_.SoundSpeed)
	return upstream_;
      else if(v<downstream_.SoundSpeed+downstream_.Velocity.x)
	return downstream_;
      else{
	const double u = (v-upstream_.SoundSpeed)*2/(1+g_);
	const double c = upstream_.SoundSpeed+(g_-1)*u/2;
	const double p = upstream_.Pressure*
	  pow(c/upstream_.SoundSpeed,2*g_/(g_-1));
	const double d = disen(p,upstream_.Density,
			       upstream_.Pressure, g_);
	return CalcPrimitive(d,p,
			     Vector2D(u,upstream_.Velocity.y),
			     IdealGas(g_));
      }
    }

  private:
    const double g_;
    const Primitive upstream_;
    const Primitive downstream_;
  };
}

/*class IsenProf: public HydroProf
  {
  public:

  IsenProf(double p2,
  Primitive const& upstream,
  double g):
  g_(g),
  upstream_(x_rest_frame(upstream)),
  downstream_
  (move_on_isentrope(p2, upstream,g)) {}

  Primitive getHydroVars(double v) const
  {
  if(v>upstream_.SoundSpeed)
  return upstream_;
  else if(v<downstream_.SoundSpeed)
  return downstream_;
  else{
  const double u = (v-upstream_.SoundSpeed)*2/(1+g_);
  const double c = upstream_.SoundSpeed+(g_-1)*u/2;
  const double p = upstream_.Pressure*
  pow(c/upstream_.SoundSpeed,2*g_/(g_-1));
  const double d = disen(p,upstream_.Density,
  upstream_.Pressure, g_);
  return CalcPrimitive(d,p,
  Vector2D(u,upstream_.Velocity.y),
  IdealGas(g_));
  }
  }

  private:
  const double g_;
  const Primitive upstream_;
  const Primitive downstream_;
  };*/

namespace{
  class GeneralProf
  {
  public:

    GeneralProf
    (double p2,
     Primitive const& upstream,
     double g):
      hugo_(max(upstream.Pressure,p2),upstream,g),
      isen_(min(upstream.Pressure,p2),upstream,g),
      hugo_flag_(p2>upstream.Pressure) {}

    Primitive getHydroVars(double v) const
    {
      if(hugo_flag_)
	return hugo_.getHydroVars(v);
      else
	return isen_.getHydroVars(v);
    }

  private:

    const ShockProf hugo_;
    const IsenProf isen_;
    const bool hugo_flag_;
  };
}

/*class GeneralProf: public HydroProf
  {
  public:

  GeneralProf(double p2,
  Primitive const& upstream,
  double g):
  hugo_flag_(p2>upstream.Pressure),
  hugo_(max(upstream.Pressure,p2),upstream,g),
  isen_(min(upstream.Pressure,p2),upstream,g) {}

  Primitive getHydroVars(double v) const
  {
  if(hugo_flag_)
  return hugo_.getHydroVars(v);
  else
  return isen_.getHydroVars(v);
  }

  private:

  const bool hugo_flag_;
  const ShockProf hugo_;
  const IsenProf isen_;
  };*/

namespace {
  class RiemannProfile
  {
  public:
    RiemannProfile(Primitive const& left,
		   Primitive const& right,
		   double g):
      contact_(calc_contact_pv(left,right,g)),
      left_prof_(contact_.pressure,
		 left, g),
      right_prof_(contact_.pressure,
		  right, g),
      left_(left), right_(right) {}

    Primitive getHydroVars(double v) const
    {
      if(v>contact_.velocity){
	Primitive res = right_prof_.getHydroVars
	  (v-right_.Velocity.x);
	res.Velocity.x += right_.Velocity.x;
	return res;
      }
      else{
	Primitive res = left_prof_.getHydroVars
	  (left_.Velocity.x - v);
	res.Velocity.x = left_.Velocity.x -
	  res.Velocity.x;
	return res;
      }
    }

  private:
    const PVContact contact_;
    const GeneralProf left_prof_;
    const GeneralProf right_prof_;
    const Primitive left_;
    const Primitive right_;
  };
}

/*class RiemannProfile//: public HydroProf
  {
  public:

  RiemannProfile(Primitive const& left,
  Primitive const& right,
  double g):
  contact_(calc_contact_pv(left,right,g)),
  left_prof_(contact_.pressure,
  left, g),
  right_prof_(contact_.pressure,
  right, g),
  left_(left), right_(right) {}

  Primitive getHydroVars(double v) const
  {
  if(v>contact_.velocity){
  Primitive res = right_prof_.getHydroVars
  (v-right_.Velocity.x);
  res.Velocity.x += right_.Velocity.x;
  return res;
  }
  else{
  Primitive res = left_prof_.getHydroVars
  (left_.Velocity.x - v);
  res.Velocity.x = left_.Velocity.x -
  res.Velocity.x;
  return res;
  }
  }

  private:
  const PVContact contact_;
  const GeneralProf left_prof_;
  const GeneralProf right_prof_;
  const Primitive left_;
  const Primitive right_;
  };*/

ERSIG::ERSIG(double g,
	     string const& vacuum_behaviour):
  g_(g), vacuum_behaviour_(vacuum_behaviour) {}

Conserved ERSIG::operator()
(Primitive const& left,
 Primitive const& right,
 double velocity) const
{
  try{
    const Primitive inter =
      RiemannProfile(left,right,g_).getHydroVars(velocity);

    const Conserved res = Primitive2Flux(inter,Vector2D(1,0)) -
      velocity*Primitive2Conserved(inter);

    return res;
  }
  catch(UniversalError& eo){
    if(eo.GetErrorMessage()=="Vacuum formed"){
      if(vacuum_behaviour_=="zero flux")
	return Conserved();
    }

    eo.AddEntry("ERSIG::Solve data starts here",0);
    eo.AddEntry("left density",left.Density);
    eo.AddEntry("left pressure",left.Pressure);
    eo.AddEntry("left x velocity",left.Velocity.x);
    eo.AddEntry("left y velocity",left.Velocity.y);
    eo.AddEntry("right density",right.Density);
    eo.AddEntry("right pressure",right.Pressure);
    eo.AddEntry("right x velocity",right.Velocity.x);
    eo.AddEntry("right y velocity",right.Velocity.y);
    throw;
  }
}
