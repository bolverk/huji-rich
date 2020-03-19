#ifndef SIMPLE_CFL_1D_HPP
#define SIMPLE_CFL_1D_HPP 1

#include "time_step_function_1d.hpp"

class SimpleCFL1D: public TimeStepFunction1D
{
public:

  explicit SimpleCFL1D(double cfl);

  double operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos) const;

private:
  double cfl_;
};

#endif // SIMPLE_CFL_1D_HPP
