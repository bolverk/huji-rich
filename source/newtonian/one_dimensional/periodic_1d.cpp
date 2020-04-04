#include "periodic_1d.hpp"
#include "../../misc/universal_error.hpp"
#include <cmath>

namespace {
  Primitive cc2primitive(const ComputationalCell& cc,
			 const EquationOfState& eos)
  {
    Primitive res;
    res.Density = cc.density;
    res.Pressure = cc.pressure;
    res.Velocity = cc.velocity;
    res.Energy = eos.dp2e(cc.density,
			  cc.pressure);
    res.SoundSpeed = eos.dp2c(cc.density,
			      cc.pressure);
    return res;
  }

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
   const size_t i) const
{
  if(i==0||i==vertex_velocity.size()-1){
    const Primitive left = cc2primitive(ss.getCells().back(), eos);
    const Primitive right = cc2primitive(ss.getCells().front(), eos);
    const double vv =
      0.5*(vertex_velocity.front() + vertex_velocity.back());
    return conserved2extensive
      (rs(left, right, vv),
       ss.getCells().back(),
       ss.getCells().front());
  }
  else
    throw UniversalError("Error in Periodic1D::CalcFlux \n Applied to bulk of grid");
}
