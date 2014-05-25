#include "profile_1d.hpp"

Profile1D::Profile1D(SpatialDistribution1D const& prof_1d):
  prof_1d_(prof_1d) {}

double Profile1D::EvalAt(Vector2D const& r) const
{
  return prof_1d_.EvalAt(r.x);
}
