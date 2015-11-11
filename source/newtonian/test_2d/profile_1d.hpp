/*! \file profile_1d.hpp
  \brief Embeds a one dimensional spatial distribution in a two dimensional simulation
  \author Almog Yalinewich
 */

#ifndef PROFILE_1D_HPP
#define PROFILE_1D_HPP 1

#include "../one_dimensional/spatial_distribution1d.hpp"
#include "../two_dimensional/spatial_distribution2d.hpp"
#include "../../tessellation/tessellation.hpp"

/*! \brief A wrapper that converts a 1d spatial profile to 2d.
  \details This distribution will only depend on the x component of the coordinate
 */
class Profile1D: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param prof_1d One dimensional profile
   */
  explicit Profile1D(SpatialDistribution1D const& prof_1d);

  double operator()(const Vector2D& r) const;

private:

  SpatialDistribution1D const& prof_1d_;
};

#endif // PROFILE_1D_HPP
