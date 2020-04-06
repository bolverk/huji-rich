#include "periodic_1d.hpp"
#include "../../misc/universal_error.hpp"
#include "flux_conversion.hpp"
#include <cmath>

Extensive Periodic1D::operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const bool /*side*/) const
{
  const Primitive left = cc2primitive(ss.getCells().back(), eos);
  const Primitive right = cc2primitive(ss.getCells().front(), eos);
  const double vv =
    0.5*(vertex_velocity.front() + vertex_velocity.back());
  return flux2extensive
    (rs(left, right, vv),
     ss.getCells().back(),
     ss.getCells().front());
}
