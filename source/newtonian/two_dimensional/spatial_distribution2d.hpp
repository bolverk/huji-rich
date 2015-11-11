/*! \file spatial_distribution2d.hpp
  \brief Spatial distribution for initial conditions
  \author Almog Yalinewich
 */

#ifndef SPATIAL_DISTRIBUTION_2D_HPP
#define SPATIAL_DISTRIBUTION_2D_HPP 1

#include "../../tessellation/tessellation.hpp"

//! \brief Spatial distribution for initial conditions
class SpatialDistribution
{
public:

  /*! \brief Evaluates the function
    \param point Position
    \return The value of the function at the point
   */
  virtual double operator()(const Vector2D& point) const = 0;

  virtual ~SpatialDistribution(void);
};

#endif // SPATIAL_DISTRIBUTION_2D_HPP
