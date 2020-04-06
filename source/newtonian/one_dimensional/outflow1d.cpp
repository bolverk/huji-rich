#include "outflow1d.hpp"
#include "../../misc/universal_error.hpp"
#include "flux_conversion.hpp"

Extensive Outflow::operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const bool side) const
{
  const vector<double> vertices = ss.getVertices();
  if(side){
    const ComputationalCell& cell = ss.getCells().back();
    const Primitive ghost = cc2primitive(cell, eos);
    const double vv = vertex_velocity[vertices.size()-1];
    return flux2extensive
      (rs(ghost,ghost,vv),
       cell);
  }
  else{
    const ComputationalCell& cell = ss.getCells().front();
    const Primitive ghost = cc2primitive(cell, eos);
    const double vv = vertex_velocity[0];
    return flux2extensive
      (rs(ghost, ghost, vv),
       cell);
  }
}
