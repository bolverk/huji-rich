/*! \file CircularRotation.hpp
  \brief Point motion that rotates the inner and outer circular points with a given angular velocity
  \author Elad Steinberg
*/

#ifndef CIRCULAR_ROTATION_HPP
#define CIRCULAR_ROTATION_HPP 1

#define _USE_MATH_DEFINES
#include "../point_motion.hpp"
#include <cmath>

//! \brief Class that calculates the angular frequency as function of radius
class OmegaFunction
{
public:
	virtual double CalcOmega(Vector2D const& point) const = 0;
};

//! \brief Class that calculates the angular frequency for a keplerian orbit
class KeplerianOmega : public OmegaFunction
{
public:
	/*! \brief Class constructor
	\param Mass The mass of the central object
	\param RigidMin The radius which below it rigid body rotation is applied
	\param RigidMax The radius which above it rigid body rotation is applied
	*/
	KeplerianOmega(double Mass, double RigidMin, double RigidMax);

	double CalcOmega(Vector2D const& point) const;
private:
	double mass_, RigidMin_, RigidMax_;
};
//! \brief Class that calculates the angular frequency for the Yee vortex
class YeeOmega : public OmegaFunction
{
public:
	/*! \brief Class constructor
	\param beta The beta parameter of the vortex
	*/
	YeeOmega(double beta) :beta_(beta){}
	double CalcOmega(Vector2D const& point)const
	{
		const double r = abs(point);
		return beta_*exp(0.5 - 0.5*r*r) / (2 * M_PI);
	}
private:
	const double beta_;
};

//! \brief Motion of mesh generating points in concentric circles, currently only works on second order time integration
class CircularRotation: public PointMotion
{
public:
	
  /*! \brief Class constructor
    \param func The function that calculates w(R)
  */
	CircularRotation(OmegaFunction const& func);


  Vector2D CalcVelocity(int index,
			Tessellation const& tess,
			vector<Primitive> const& cells,
			double time);
 
  void ApplyFix(Tessellation const& tess, vector<Primitive> const& cells, double time,
	  vector<CustomEvolution*> &cevolve, const vector<vector<double> >& tracers, double dt, vector < Vector2D >
	  & velocities);

private:
  bool evencall_;
  OmegaFunction const& omega_;
};




#endif // CIRCULAR_ROTATION_HPP
