#include "eulerian1d.hpp"

double Eulerian1D::operator()
(size_t /*i*/, vector<double> const& /*vp*/,
 vector<ComputationalCell> const& /*hv*/) const
{
  return 0;
}
