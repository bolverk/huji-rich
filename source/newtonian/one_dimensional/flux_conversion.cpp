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
