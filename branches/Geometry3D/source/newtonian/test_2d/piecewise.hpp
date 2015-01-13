/*! \file piecewise.hpp
  \author Almog Yalinewich
  \brief Piecewise spatial distribution
 */

#ifndef PIECEWISE_HPP
#define PIECEWISE_HPP 1

#include "../two_dimensional/spatial_distribution2d.hpp"
#include "../../tessellation/shape_2d.hpp"

/*! \brief A piecewise distribution function
  \details Uses one distribution function if a point is withing a specified shape, and another if it is outside
 */
class Piecewise: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param shape A shape
    \param inside Function to use when inside the shape
    \param outside Function to use when outside the shape
   */
  Piecewise(const Shape2D& shape,
	    const SpatialDistribution& inside,
	    const SpatialDistribution& outside);

  double operator()(const Vector2D& point) const;

private:
  const Shape2D& shape_;
  const SpatialDistribution& inside_;
  const SpatialDistribution& outside_;
};

#endif // PIECEWISE_HPP
