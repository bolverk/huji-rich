#include <cassert>
#include "hdsim_3d.hpp"
#include "extensive_generator.hpp"

HDSim3D::HDSim3D(const Tessellation3D& tess,
		 const vector<ComputationalCell>& cells,
		 const EquationOfState& eos,
		 const PointMotion3D& pm,
		 const TimeStepCalculator& tsc,
		 const FluxCalculator& fc,
		 const CellUpdater& cu):
  tess_(tess), eos_(eos), cells_(cells),
  extensive_(serial_generate(ExtensiveGenerator(cells,tess,eos))),
  pm_(pm), tsc_(tsc), fc_(fc), cu_(cu)
{
  assert(tess.GetPointNo()==cells.size());
}
