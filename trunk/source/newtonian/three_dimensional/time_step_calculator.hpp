#ifndef TIME_STEP_CALCULATOR_HPP
#define TIME_STEP_CALCULATOR_HPP 1

class TimeStepCalculator
{
public:

  virtual double operator()
  (const Tessellation3D& tess,
   const vector<ComputationalCell>& cells) const = 0;
};

#endif  // TIME_STEP_CALCULATOR_HPP
