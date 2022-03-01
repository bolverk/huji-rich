#include "equation_of_state.hpp"

EquationOfState::~EquationOfState(void) {}

double EquationOfState::dT2cv(double const d, double const T,
	  tvector const& tracers, vector<string> const& tracernames) const
{
    throw UniversalError("dT2cv not implemented");
}

double EquationOfState::de2T(double const d, double const e,
	  tvector const& tracers, vector<string> const& tracernames) const
{
    throw UniversalError("de2T not implemented");
}