#include "equation_of_state.hpp"
#include "source/misc/universal_error.hpp"

EquationOfState::~EquationOfState(void) {}

double EquationOfState::dT2cv(double const /*d*/, double const /*T*/,
	  tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
    throw UniversalError("dT2cv not implemented");
    return 0;
}

double EquationOfState::de2T(double const /*d*/, double const /*e*/,
	  tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
    throw UniversalError("de2T not implemented");
    return 0;
}

double EquationOfState::dT2e(double const /*d*/, double const /*T*/,
	  tvector const& /*tracers*/, vector<string> const& /*tracernames*/) const
{
   throw UniversalError("dT2e not implemented");
   return 0;
}