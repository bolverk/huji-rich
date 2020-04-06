#include "periodic_1d.hpp"
#include "../../misc/universal_error.hpp"
#include "flux_conversion.hpp"
#include <cmath>

namespace {
  Extensive conserved2extensive
  (const Conserved& flux,
   const ComputationalCell& left,
   const ComputationalCell& right)
  {
    Extensive res;
    res.mass = flux.Mass;
    res.momentum = flux.Momentum;
    res.energy = flux.Energy;
    for(size_t i=0;i<left.tracers.size();++i)
      res.tracers.push_back
	(res.mass*(res.mass > 0 ?
		   left.tracers.at(i) :
		   right.tracers.at(i)));
    return res;
  }
}

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
  return conserved2extensive
    (rs(left, right, vv),
     ss.getCells().back(),
     ss.getCells().front());
}
