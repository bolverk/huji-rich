#include "ResetDump.hpp"

ResetDump::ResetDump(void):snapshot(vector<Vector2D> (),vector<Primitive> ()),
	tracers(vector<vector<double> > ()),time(0),cycle(0),coldflows(false),densityfloor(false),a(0),b(0),
	densitymin(0),pressuremin(0),procmesh(vector<Vector2D>()),cevolve(vector<size_t> ())
{}

ResetDump::~ResetDump(void)
{
clear();
}

void ResetDump::clear(void)
{
	tracers.clear();
	snapshot.cells.clear();
	snapshot.mesh_points.clear();
}

HydroSnapshot::HydroSnapshot
(vector<Vector2D> const& mesh_points_i,
 vector<Primitive> const& cells_i):
  mesh_points(mesh_points_i),
  cells(cells_i) {}

HydroSnapshot::HydroSnapshot(void):
  mesh_points(vector<Vector2D>()),
  cells(vector<Primitive>()) {}
