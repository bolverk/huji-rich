/*! \file Circle2D.hpp
  \brief A spatial distribution with a Circle step function
  \author Elad Steinberg
 */

#ifndef CIRCLE2D_HPP
#define CIRCLE2D_HPP 1

#include "../spatial_distribution2d.hpp"

//! \brief A spatial distribution with a square step function
class Circle2D: public SpatialDistribution
{
private:

  double _R,_in,_out;
  Vector2D center;

public:

  /*! \brief Class constructor
    \param xc X coordinate of the circle center
    \param yc Y coordinate of the circle center
    \param R Radius of the circle
    \param in Value inside the circle
    \param out Value outside the circle
   */
  Circle2D(double xc,double yc,double R, double in,double out);
  ~Circle2D();
  double EvalAt(Vector2D const& r) const;
};

#endif // CIRCLE2D_HPP
