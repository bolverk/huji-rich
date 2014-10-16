#include <cmath>
#include "ideal_gas.hpp"
#include "../../misc/universal_error.hpp"

IdealGas::IdealGas(double AdiabaticIndex):
  _g(AdiabaticIndex) {}

double IdealGas::getAdiabaticIndex(void) const
{
  return _g;
}

double IdealGas::dp2e(double d, double p) const
{
  return p/d/(_g-1);
}

double IdealGas::de2p(double d, double e) const
{
  return (_g-1)*e*d;
}

namespace {
  UniversalError imaginary_speed_of_sound(double g,
					  double p,
					  double d)
  {
    UniversalError res("Speed of sound came out imaginary");
    res.AddEntry("adiabatic index",g);
    res.AddEntry("density",d);
    res.AddEntry("pressure",p);
    return res;
  }
}

double IdealGas::dp2c(double d, double p) const
{
  if(_g<=0||p<=0||d<=0||is_nan(d)||is_nan(p))
    throw imaginary_speed_of_sound(_g,p,d);
  return sqrt(_g*p/d);
}

double IdealGas::de2c(double d, double e) const
{
  double p = de2p(d, e);
  return sqrt(_g*p/d);
}

double IdealGas::dp2s(double d, double p) const
{
	return p*pow(d,-_g);
}

double IdealGas::sd2p(double s, double d) const
{
	return s*pow(d,_g);
}
