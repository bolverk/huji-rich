#include "eulerian1d.hpp"

double Eulerian1D::operator()
(int /*i*/, vector<double> const& /*vp*/,
 vector<Primitive> const& /*hv*/) const
{
  return 0;
}
