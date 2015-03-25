#include "profile_1d.hpp"

Profile1D::Profile1D(SpatialDistribution1D const& prof_1d):
  prof_1d_(prof_1d) {}

double Profile1D::operator()(Vector2D const& r) const
{
  return prof_1d_(r.x);
}
