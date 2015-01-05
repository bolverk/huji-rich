#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include "conserved_3d.hpp"

class FluxCalculator
{
public:

  virtual vector<Conserved3D> operator()
  (const Tessellation3D& tess,
   const vector<ComputationalCell>& cells) = 0;
};

#endif // FLUX_CALCULATOR_HPP
