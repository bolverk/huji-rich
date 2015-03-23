#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include <vector>
#include "extensive.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"

using std::vector;

class FluxCalculator
{
public:

  virtual vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& face_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const double time,
   const double dt) const = 0;

  virtual ~FluxCalculator(void);
};

#endif // FLUX_CALCULATOR_HPP
