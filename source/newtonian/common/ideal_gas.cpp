#include <cmath>
#include "ideal_gas.hpp"
#include "../../misc/universal_error.hpp"
#include <limits>

IdealGas::IdealGas(double const AdiabaticIndex) : g_(AdiabaticIndex), f_(std::numeric_limits<double>::signaling_NaN()), 
	beta_(std::numeric_limits<double>::signaling_NaN()), mu_(std::numeric_limits<double>::signaling_NaN()) {}

IdealGas::IdealGas(double const AdiabaticIndex, double const f, double const beta, double const mu) :
	g_(AdiabaticIndex), f_(f), beta_(beta), mu_(mu) {}

double IdealGas::getAdiabaticIndex(void) const
{
	return g_;
}

double IdealGas::dp2e(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
#ifdef RICH_DEBUG
	if (!(g_ > 0) || !(p > 0) || !(d > 0))
	{
		UniversalError eo("Negative quantity in ideal gas dp2e");
		eo.addEntry("Density", d);
		eo.addEntry("Pressure", p);
		eo.addEntry("Gamma index", g_);
		throw eo;
	}
#endif //RICH_DEBUG

	return p / d / (g_ - 1);
}

double IdealGas::de2p(double d, double e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
#ifdef RICH_DEBUG
	if (!(g_ > 0) || !(e > 0) || !(d > 0))
	{
		UniversalError eo("Negative quantity in ideal gas de2p");
		eo.addEntry("Density", d);
		eo.addEntry("Energy", e);
		eo.addEntry("Gamma index", g_);
		throw eo;
	}
#endif // RICH_DEBUG

	if (e < 0)
		throw UniversalError("Negative thermal energy");
	return (g_ - 1)*e*d;
}

namespace {
	/*
	UniversalError imaginary_speed_of_sound(double g,
						double p,
						double d)
	{
	  UniversalError res("Speed of sound came out imaginary");
	  res.addEntry("adiabatic index",g);
	  res.addEntry("density",d);
	  res.addEntry("pressure",p);
	  return res;
	}
	*/
}

double IdealGas::dp2c(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
#ifdef RICH_DEBUG
	if (!(g_ > 0) || !(p > 0) || !(d > 0))
	{
		UniversalError eo("Negative quantity in ideal gas");
		eo.addEntry("Density", d);
		eo.addEntry("Pressure", p);
		eo.addEntry("Gamma index", g_);
		throw eo;
	}
#endif //RICH_DEBUG
	assert(!(p < 0));
	assert(!(d < 0));
	assert(!(g_ < 0));
	return sqrt(g_*p / d);
}

double IdealGas::de2c(double d, double e, tvector const& tracers, vector<string> const& tracernames) const
{
	double p = de2p(d, e, tracers, tracernames);
	return sqrt(g_*p / d);
}

double IdealGas::dp2s(double d, double p, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return p * pow(d, -g_);
}

double IdealGas::sd2p(double s, double d, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	if ((!(s > 0) || (!(d > 0))))
	{
		UniversalError eo("Negative quantaties in entropy calc");
		eo.addEntry("entropy", s);
		eo.addEntry("density", d);
		throw eo;
	}
	assert(s > 0 && d > 0);
	return s * pow(d, g_);
}

double IdealGas::dT2cv(double const d, double const T, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return f_ * beta_ * std::pow(T, beta_ - 1) * std::pow(d, -mu_);
}

double IdealGas::de2T(double const d, double const e, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return std::pow(std::pow(d, mu_) * e / f_, 1.0 / beta_);
}

double IdealGas::dT2e(double const d, double const T, tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
	return f_ * std::pow(T, beta_) * std::pow(d, -mu_);
}