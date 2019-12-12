#include <cmath>
#include "ideal_gas.hpp"
#include "../../misc/universal_error.hpp"

IdealGas::IdealGas(double AdiabaticIndex):
  g_(AdiabaticIndex) {}

double IdealGas::getAdiabaticIndex(void) const
{
  return g_;
}

double IdealGas::dp2e(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
  return p/d/(g_-1);
}

double IdealGas::de2p(double d, double e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if (e < 0)
		throw UniversalError("Negative thermal energy");
  return (g_-1)*e*d;
}

namespace {
  /*
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
  */
}

double IdealGas::dp2c(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if (!(d > 0) || !(p > 0))
	{
		UniversalError eo("Bad cs");
		eo.AddEntry("pressure", p);
		eo.AddEntry("density", d);
		throw eo;
	}
  assert(g_>0 && p>0 && d>0);
  return sqrt(g_*p/d);
}

double IdealGas::de2c(double d, double e, tvector const& tracers, vector<string> const& tracernames) const
{
  double p = de2p(d, e, tracers,tracernames);
  return sqrt(g_*p/d);
}

double IdealGas::dp2s(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
  return p*pow(d,-g_);
}

double IdealGas::sd2p(double s, double d, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
  assert(s>0 && d>0);
  return s*pow(d,g_);
}
