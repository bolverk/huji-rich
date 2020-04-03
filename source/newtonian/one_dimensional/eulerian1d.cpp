#include "eulerian1d.hpp"

double Eulerian1D::operator()
(int /*i*/, vector<double> const& /*vp*/,
 vector<ComputationalCell> const& /*hv*/) const
{
  return 0;
}
