#include "simple_cfl.hpp"
#include "hydrodynamics_2d.hpp"

SimpleCFL::SimpleCFL(const double cfl): cfl_(cfl) {}

double SimpleCFL::operator()(const Tessellation& tess,
			     const vector<Primitive>& cells,
			     const vector<Vector2D>& point_velocities,
			     const HydroBoundaryConditions& hbc,
			     const double time,
			     const vector<CustomEvolution*>& custom_evolution)
{
  return cfl_*CalcTimeStep(tess,cells,
			   tess.calc_edge_velocities(&hbc,point_velocities,time),
			   hbc,time,custom_evolution);
}
