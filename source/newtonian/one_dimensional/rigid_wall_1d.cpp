#include "rigid_wall_1d.hpp"
#include "../../misc/universal_error.hpp"
#include "flux_conversion.hpp"

namespace {
  Primitive reverse_velocity(const Primitive& p)
  {
    Primitive res = p;
    res.Velocity.x *= -1;
    return res;
  }
}

Extensive RigidWall1D::operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs, 
   const vector<double>& vertex_velocity,
   const bool side) const
{
  const vector<double>& vertices = ss.getVertices();
  if(side){
    const ComputationalCell& cell = ss.getCells().back();
    const Primitive left = cc2primitive(cell, eos);
    const Primitive right = reverse_velocity(left);
    const double vv = vertex_velocity.at(vertices.size()-1);
    return flux2extensive
      (rs(left, right, vv),
       cell);
  }
  else{
    const ComputationalCell& cell = ss.getCells().front();
    const Primitive right = cc2primitive(cell, eos);
    const Primitive left = reverse_velocity(right);
    const double vv = vertex_velocity[0];
    return flux2extensive
      (rs(left, right, vv),
       cell);
  }
}
