#ifndef TIME_STEP_CALCULATOR_HPP
#define TIME_STEP_CALCULATOR_HPP 1

#include "../../3D/Tessellation/Tessellation3D.hpp"
#include "computational_cell.hpp"

class TimeStepCalculator
{
public:

  virtual double operator()
  (const Tessellation3D& tess,
   const vector<ComputationalCell>& cells) const = 0;

  virtual ~TimeStepCalculator(void);
};

#endif  // TIME_STEP_CALCULATOR_HPP
