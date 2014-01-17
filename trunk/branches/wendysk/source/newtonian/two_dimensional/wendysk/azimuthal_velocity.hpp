#ifndef AZIMUTHAL_VELOCITY_HPP
#define AZIMUTHAL_VELOCITY_HPP 1

#include "../spatial_distribution2d.hpp"
#include "../../../misc/func_1_var.hpp"

class AzimuthalVelocity: public SpatialDistribution
{
public:

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
