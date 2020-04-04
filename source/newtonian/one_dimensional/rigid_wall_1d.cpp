#include "rigid_wall_1d.hpp"
#include "../../misc/universal_error.hpp"

namespace {
  Primitive reverse_velocity(const Primitive& p)
  {
    Primitive res = p;
    res.Velocity.x *= -1;
    return res;
  }

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
  (const Conserved& conserved,
   const ComputationalCell& cell)
  {
    Extensive res;
    res.mass = conserved.Mass;
    res.momentum = conserved.Momentum;
    res.energy = conserved.Energy;
    res.tracers = vector<double>(cell.tracers.size(),0);
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
    return conserved2extensive
      (rs(left, right, vv),
       cell);
  }
  else{
    const ComputationalCell& cell = ss.getCells().front();
    const Primitive right = cc2primitive(cell, eos);
    const Primitive left = reverse_velocity(right);
    const double vv = vertex_velocity[0];
    return conserved2extensive
      (rs(left, right, vv),
       cell);
  }
}
