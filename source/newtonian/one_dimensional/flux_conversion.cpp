#include "flux_conversion.hpp"

Primitive cc2primitive
(const ComputationalCell& cell,
 const EquationOfState& eos)
{
  Primitive res;
  res.Density = cell.density;
  res.Pressure = cell.pressure;
  res.Velocity = cell.velocity;
  res.Energy = eos.dp2e(cell.density,
			cell.pressure);
  res.SoundSpeed = eos.dp2c(cell.density,
			    cell.pressure);
  return res;
}

vector<Primitive> ccs2primitives
(const vector<ComputationalCell>& cells,
 const EquationOfState& eos)
{
  vector<Primitive> res(cells.size());
  for(size_t i=0;i<res.size();++i)
    res.at(i) = cc2primitive(cells.at(i), eos);
  return res;
}

Extensive flux2extensive
(const Conserved& flux,
 const ComputationalCell& donor)
{
  Extensive res;
  res.mass = flux.Mass;
  res.momentum = flux.Momentum;
  res.energy = flux.Energy;
  for(size_t i=0;i<donor.tracers.size();++i)
    res.tracers.push_back(res.mass*donor.tracers.at(i));
  return res;
}

Extensive flux2extensive
(const Conserved& flux,
 const ComputationalCell& left,
 const ComputationalCell& right)
{
  return flux2extensive
    (flux,
     flux.Mass > 0 ? left : right);
}
