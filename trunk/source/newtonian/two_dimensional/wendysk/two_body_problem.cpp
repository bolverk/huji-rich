#include <cmath>
#include "two_body_problem.hpp"
#include "../../../misc/bisection.hpp"

TwoBodyProblem::GeometricalPoint::GeometricalPoint
(Vector2D const& position_i,
 Vector2D const& velocity_i):
  position(position_i),
  velocity(velocity_i) {}

TwoBodyProblem::GeometricalPoint::~GeometricalPoint(void) {}

TwoBodyProblem::Particle::Particle(Vector2D const& position_i,
				   Vector2D const& velocity_i,
				   double mass_i):
  GeometricalPoint(position_i, velocity_i),
  mass(mass_i) {}

TwoBodyProblem::Particle::Particle(Particle const& source):
  GeometricalPoint(source.position,
		   source.velocity),
  mass(source.mass) {}

namespace {

  class KeplerEquation: public Func1Var
  {
  public:

    KeplerEquation(double eccentricity,
		   double rhs):
      eccentricity_(eccentricity),
      rhs_(rhs) {}

    double eval(double x) const
    {
      return x-eccentricity_*sin(x)-rhs_;
    }

  private:
    const double eccentricity_;
    const double rhs_;
  };

  double solve_kepler_equation(double eccentricity,
			       double rhs)
  {
    return bisection(KeplerEquation(eccentricity,
				    rhs),
		     rhs-eccentricity,
		     rhs+eccentricity);
  }

  double cross_z_comp(Vector2D const& v1,
		      Vector2D const& v2)
  {
    return v1.x*v2.y-v1.y*v2.x;
  }

  Vector2D mass_average(double m1,
			Vector2D v1,
			double m2,
			Vector2D v2)
  {
    return (m1*v1+m2*v2)/(m1+m2);
  }

  TwoBodyProblem::GeometricalPoint calc_center_of_mass
  (pair<TwoBodyProblem::Particle,TwoBodyProblem::Particle> const& pp)
  {
    return TwoBodyProblem::GeometricalPoint
      (mass_average(pp.first.mass,
		    pp.first.position,
		    pp.second.mass,
		    pp.second.position),
       mass_average(pp.first.mass,
		    pp.first.velocity,
		    pp.second.mass,
		    pp.second.velocity));
  }

  TwoBodyProblem::GeometricalPoint calc_deviator
  (pair<TwoBodyProblem::Particle,TwoBodyProblem::Particle> const& pp)
  {
    return TwoBodyProblem::GeometricalPoint(pp.second.position-
					    pp.first.position,
					        pp.second.velocity-
					    pp.first.velocity);
  }

  double calc_eccentricity(double e,
			   double l,
			   double k)
  {
    return sqrt(1+2*e*pow(l,2)/(pow(k,2)));
  }

  double calc_eccentric_anomaly(double angle,
				double eccentricity)
  {
    return acos((eccentricity+cos(angle))/
		(1+eccentricity*cos(angle)));
  }

  double eccentric_anomaly_to_angle(double eccentricity,
				    double eccentric_anomaly)
  {
    const double tane2 = tan(0.5*eccentric_anomaly);
    const double eps = eccentricity;
    return 2.*atan(sqrt(((1.+eps)/(1.-eps)))*tane2);
  }

  double calc_energy(TwoBodyProblem::GeometricalPoint const& deviator,
		     double k)
  {
    const double kinetic = 0.5*pow(abs(deviator.velocity),2);
    const double potential = -k/abs(deviator.position);
    return kinetic+potential;
  }

  double calc_angular_momentum(TwoBodyProblem::GeometricalPoint const& deviator)
  {
    return cross_z_comp(deviator.position,
			deviator.velocity);
  }

  Vector2D calc_lrl(TwoBodyProblem::GeometricalPoint const& deviator,
		    double k)
  {
    return deviator.position*ScalarProd(deviator.velocity,
					deviator.velocity)-
      deviator.velocity*ScalarProd(deviator.position,
				   deviator.velocity)-
      deviator.position*k/abs(deviator.position);
  }

  double calc_time_scale(double angular_momentum,
			 double k,
			 double eccentricity)
  {
    return (pow(angular_momentum,3)/pow(k,2))/
      pow(1.-pow(eccentricity,2),1.5);
  }

  double calc_initial_phase(TwoBodyProblem::GeometricalPoint const& deviator,
			    Vector2D const& lrl,
			    double eccentricity,
			    double t2timescale)
  {
    const double initial_angle = CalcAngle(deviator.position,
					   lrl);
    const double initial_eccentric_anomaly =
      calc_eccentric_anomaly(initial_angle,
			     eccentricity);
    return t2timescale-
      initial_eccentric_anomaly+
      eccentricity*sin(initial_eccentric_anomaly);
  }
}

TwoBodyProblem::TwoBodyProblem
(pair<TwoBodyProblem::Particle,TwoBodyProblem::Particle> const& pp,
 double initial_time,
 double g):
  pp_(pp.first, pp.second),
  k_(g*(pp.first.mass+pp.second.mass)),
  center_of_mass_(calc_center_of_mass(pp)),
  deviator_(calc_deviator(pp)),
  energy_(calc_energy(deviator_,k_)),
  angular_momentum_(calc_angular_momentum(deviator_)),
  lrl_(calc_lrl(deviator_,k_)),
  eccentricity_(calc_eccentricity(energy_,
				  angular_momentum_,
				  k_)),
  time_scale_(calc_time_scale(angular_momentum_,
			      k_, eccentricity_)),
  initial_phase_(calc_initial_phase(deviator_,
				    lrl_,
				    eccentricity_,
				    initial_time/time_scale_)),
  angle_offset_(acos(lrl_.x/abs(lrl_))) {}
  
pair<TwoBodyProblem::GeometricalPoint, TwoBodyProblem::GeometricalPoint> 
TwoBodyProblem::timeEvolve(double t) const
{
  const double rhs = t/time_scale_-initial_phase_;
  const double eccentric_anomaly = solve_kepler_equation
    (eccentricity_,rhs);
  const double angle = eccentric_anomaly_to_angle
    (eccentricity_, eccentric_anomaly);
  const double radius = (pow(angular_momentum_,2)/k_)/
    (1.+eccentricity_*cos(angle));
  const Vector2D deviator_pos = pol2cart(radius,
					 angle+angle_offset_);
  const Vector2D cm_pos = center_of_mass_.position+
    center_of_mass_.velocity*t;
  const Vector2D pos1 = cm_pos-deviator_pos*
    (pp_.second.mass/(pp_.first.mass+pp_.second.mass));
  const Vector2D pos2 = cm_pos+deviator_pos*
    (pp_.first.mass/(pp_.first.mass+pp_.second.mass));
  const double tangential_velocity = angular_momentum_/abs(deviator_pos);
  const double radial_velocity = 2*energy_+ 2*k_/abs(deviator_pos);
  const Vector2D deviator_vel = pol2cart(radial_velocity,
					 tangential_velocity);
  const Vector2D vel1 = center_of_mass_.velocity - deviator_vel*
    (pp_.second.mass/(pp_.first.mass+pp_.second.mass));
  const Vector2D vel2 = center_of_mass_.velocity + deviator_vel*
    (pp_.first.mass/(pp_.first.mass+pp_.second.mass));
  return pair<GeometricalPoint,GeometricalPoint>
    (GeometricalPoint(pos1,vel1),
     GeometricalPoint(pos2,vel2));
}

double TwoBodyProblem::getGravitationConstant(void) const
{
  return k_/(pp_.first.mass+pp_.second.mass);
}

pair<TwoBodyProblem::Particle,TwoBodyProblem::Particle> const& TwoBodyProblem::getInitialParticles(void) const
{
  return pp_;
}

TwoBodyProblem::~TwoBodyProblem(void) {}
