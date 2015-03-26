#include "extensive_generator.hpp"

Conserved3D calc_intensive(const ComputationalCell& cell,
			   const EquationOfState& eos)
{
  return Conserved3D(cell.density,
		     cell.density*cell.velocity,
		     cell.density*
		     (0.5*ScalarProd(cell.velocity,cell.velocity)+
		      eos.dp2e(cell.density,cell.pressure)),
		     cell.tracers);
}

ExtensiveGenerator::ExtensiveGenerator
(const vector<ComputationalCell>& cells,
 const Tessellation3D& tess,
 const EquationOfState& eos):
  cells_(cells), tess_(tess), eos_(eos) {}

size_t ExtensiveGenerator::size(void) const
{
  return cells_.size();
}

Conserved3D ExtensiveGenerator::operator[](size_t i) const
{
  return tess_.GetVolume(i)*calc_intensive(cells_[i],eos_);
}
