#include "outflow1d.hpp"
#include "../../misc/universal_error.hpp"

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
}

Extensive Outflow::operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const size_t i) const
{
  const vector<double> vertices = ss.getVertices();
  if(i==0){
    const ComputationalCell& cell = ss.getCells().front();
    const Primitive ghost = cc2primitive(cell, eos);
    const double vv = vertex_velocity[0];
    const Conserved flux = rs(ghost,ghost,vv);
    Extensive res;
    res.mass = flux.Mass;
    res.momentum = flux.Momentum;
    res.energy = flux.Energy;
    for(size_t j=0;j<cell.tracers.size();++j)
      res.tracers.push_back(res.mass*cell.tracers.at(j));
    return res;	 
  }
  else if(i==vertices.size()-1){
    const ComputationalCell& cell = ss.getCells().back();
    const Primitive ghost = cc2primitive(cell, eos);
    const double vv = vertex_velocity[vertices.size()-1];
    const Conserved flux = rs(ghost,ghost,vv);
    Extensive res;
    res.mass = flux.Mass;
    res.momentum = flux.Momentum;
    res.energy = flux.Energy;
    for(size_t j=0;j<cell.tracers.size();++j)
      res.tracers.push_back(res.mass*cell.tracers.at(j));
    return res;
  }
  else
    throw UniversalError("Index inside bulk of grid");
}
