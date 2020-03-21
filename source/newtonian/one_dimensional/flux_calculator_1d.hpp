#ifndef FLUX_CALCULATOR_1D_HPP
#define FLUX_CALCULATOR_1D_HPP 1

#include <vector>
#include "../two_dimensional/extensive.hpp"
#include "simulation_state_1d.hpp"
#include "../common/equation_of_state.hpp"

using std::vector;

class FluxCalculator1D
{
public:

  virtual vector<Extensive> operator()
  (const SimulationState1D& ss,
   const vector<double>& vertex_velocity,
   const EquationOfState& eos,
   const double dt) const = 0;

  virtual ~FluxCalculator1D(void);
};

#endif // FLUX_CALCULATOR_1D_HPP
