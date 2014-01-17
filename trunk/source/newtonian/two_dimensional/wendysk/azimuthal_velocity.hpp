/*! \file azimuthal_velocity.hpp
  \brief Decomposition of azimuthal velocity in polar coordinates to x and y components
  \author Almog Yalinewich
 */

#ifndef AZIMUTHAL_VELOCITY_HPP
#define AZIMUTHAL_VELOCITY_HPP 1

#include "../spatial_distribution2d.hpp"
#include "../../../misc/func_1_var.hpp"

//! \brief Decomposition of azimuthal velocity in polar coordinates to x and y components
class AzimuthalVelocity: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param radial Azimuthal velocity as a function of radius
    \param center Center
    \param comp either x or y component
   */
  AzimuthalVelocity(Func1Var const& radial,
		    Vector2D const& center,
		    char comp);

  double EvalAt(Vector2D const& p) const;

private:
  Func1Var const& radial_;
  const Vector2D center_;
  const char comp_;
};

#endif // AZIMUTHAL_VELOCITY_HPP
