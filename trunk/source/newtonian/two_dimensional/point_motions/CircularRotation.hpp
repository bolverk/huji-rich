/*! \file CircularRotation.hpp
  \brief Point motion that rotates the inner and outer circular points with a given angular velocity
  \author Elad Steinberg
*/

#ifndef CIRCULAR_ROTATION_HPP
#define CIRCULAR_ROTATION_HPP 1

#include "../point_motion.hpp"
#include <cmath>

//! \brief Motion of mesh generating points in concentric circles
class CircularRotation: public PointMotion
{
public:
  /*! \brief Class constructor
    \param naive The original point motion
    \param Rinner The radius which beneath it points are rotated with w_inner
    \param Router The radius which above it points are rotated with w_inner
    \param w_inner The angular velocity of the inner radius
    \param w_outer The angular velocity of the outer radius
    \param Ninner The number of points in the inner circle. The point ordering is assumed to be inner_circles,outermost_circles,other_points
    \param Nouter The number of points in the outer circle.
    \param t The time of the simulation
    \param center The center of the circles
  */
  CircularRotation(PointMotion& naive,double Rinner,double Router,double w_inner,
		   double w_outer,int Ninner,int Nouter,double t=0,
		   Vector2D const& center=Vector2D(0,0));

  Vector2D CalcVelocity(int index,
			Tessellation const& tess,
			vector<Primitive> const& cells,
			double time);

  vector<Vector2D> calcAllVelocities(Tessellation const& tess,
				     vector<Primitive> const& cells,
				     double time,vector<CustomEvolution*> & cevolve);

private:
  PointMotion& naive_;
  const double inner_radius_;
  const double outer_radius_;
  const double w_inner_;
  const double w_outer_;
  const int Ninner_;
  const int Nouter_;
  const double t0_;
  const Vector2D center_;
  vector<double> init_angles_;
  vector<double> init_R_;
  bool first_time_;
};

#endif // CIRCULAR_ROTATION_HPP
