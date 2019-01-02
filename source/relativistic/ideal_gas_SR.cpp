#include <cmath>
#include "ideal_gas_SR.hpp"
#include "../misc/universal_error.hpp"

IdealGas_SR::IdealGas_SR(double AdiabaticIndex) :
	g_(AdiabaticIndex) {}

double IdealGas_SR::getAdiabaticIndex(void) const
{
	return g_;
}

double IdealGas_SR::dp2e(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return p * g_ / (d*(g_ - 1));
}

double IdealGas_SR::de2p(double d, double e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if (e < 0)
		throw UniversalError("Negative thermal energy");
	return e * (g_ - 1)*d / g_;
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

double IdealGas_SR::dp2c(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
#ifdef RICH_DEBUG
	if (!(g_ > 0) || !(p > 0) || !(d > 0))
	{
		UniversalError eo("Negative quantity in ideal gas SR");
		eo.AddEntry("Density", d);
		eo.AddEntry("Pressure", p);
		eo.AddEntry("Gamma index", g_);
		throw eo;
	}
#endif
	assert(!(p < 0));
	assert(!(d < 0));
	assert(!(g_ < 0));
	return std::sqrt(g_*p / (d + p * g_ / ((g_ - 1))));
}

double IdealGas_SR::de2c(double d, double e, tvector const& tracers, vector<string> const& tracernames) const
{
	double p = de2p(d, e, tracers, tracernames);
	return std::sqrt(g_*p / (d*(1 + e)));
}

double IdealGas_SR::dp2s(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return p * pow(d, -g_);
}

double IdealGas_SR::sd2p(double s, double d, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	assert(s > 0 && d > 0);
	return s * pow(d, g_);
}