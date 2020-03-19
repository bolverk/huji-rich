#ifndef TIME_STEP_FUNCTION_1D
#define TIME_STEP_FUNCTION_1D 1

#include "simulation_state_1d.hpp"
#include "../common/equation_of_state.hpp"

class TimeStepFunction1D
{
public:

  virtual double operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos) const = 0;

  virtual ~TimeStepFunction1D(void);
};

#endif // TIME_STEP_FUNCTION_1D
