/*! \file Line2D.hpp
  \brief A spatial distribution with a Line step function
  \author Elad Steinberg
 */

#ifndef LINE2D_HPP
#define LINE2D_HPP 1

#include "../spatial_distribution2d.hpp"

//! \brief A spatial distribution with a square step function
class Line2D: public SpatialDistribution
{
private:

  double _a,_b,_above,_under;

public:

  /*! \brief Class constructor
    \param a Slope of the line a*x+b=y
    \param b Where the line intersects zero
    \param above Value above the line
    \param under Value under the line
   */
  Line2D(double a,double b,double above,double under);
  ~Line2D();
  double operator()(const Vector2D& r) const;
};

#endif // LINE2D_HPP
