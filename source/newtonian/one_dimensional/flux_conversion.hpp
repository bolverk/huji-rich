#ifndef FLUX_CONVERSION_HPP
#define FLUX_CONVERSION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"
#include "../common/equation_of_state.hpp"

Primitive cc2primitive
(const ComputationalCell& cell,
 const EquationOfState& eos);

vector<Primitive> ccs2primitives
(const vector<ComputationalCell>& cells,
 const EquationOfState& eos);

#endif // FLUX_CONVERSION_HPP
