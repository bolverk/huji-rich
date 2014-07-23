#include "FreeFlow.hpp"
#include "../hydrodynamics_2d.hpp"

double FreeFlow::GetMassFlux(void) const
{
	return mass_flux;
}

FreeFlow::FreeFlow(RiemannSolver const& rs,bool mass_count):
  rs_(rs),mass_count_(mass_count), mass_flux(0) {}

FreeFlow::~FreeFlow(void) {}

Primitive FreeFlow::GetBoundaryPrimitive
	(Edge const& edge,
	Tessellation const& /*Data*/,
	vector<Primitive> const& cells,double /*time*/)const
{
	if(edge.neighbors.first==-1)
		return cells[edge.neighbors.second];
	else
		return cells[edge.neighbors.first];
}

vector<double> FreeFlow::GetBoundaryTracers
	(Edge const& edge,
	Tessellation const& /*Data*/,
	 vector<vector<double> > const& tracers,double
	/*time*/)const
{
	if(edge.neighbors.first==-1)
		return tracers[edge.neighbors.second];
	else
		return tracers[edge.neighbors.first];
}

namespace {
int calc_ci(Edge const& edge)
{
	if((edge.neighbors.first==-1)&&
		(edge.neighbors.second!=-1))
		return 1;
	else if((edge.neighbors.second==-1)&&
		edge.neighbors.first!=-1)
		return 0;
	else
	  throw UniversalError("Boundary condition called for bulk cell");
}
}

Conserved FreeFlow::CalcFlux
(Tessellation const& tessellation,
 vector<Primitive> const& cells,Vector2D const& edge_velocity,
 Edge const& edge,
 SpatialReconstruction const& /*interp*/,double /*dt*/,
 double /*time*/) const
{
	const int ci = calc_ci(edge);
	const Vector2D p = Parallel(edge);
	const Vector2D n = Normal(edge, tessellation);
	const Primitive ghost=cells[pair_member(edge.neighbors,ci)];
	vector<Primitive> states(2);
	for(int i=0;i<2;i++){
	  states[i] = ghost;
		states[i].Velocity.Set
			(Projection(states[i].Velocity, n),
			Projection(states[i].Velocity, p));
	}
	Conserved res = rs_.Solve(states[0], states[1],Projection(edge_velocity,n));
	res.Momentum = res.Momentum.x*n/abs(n) +
		res.Momentum.y*p/abs(p);
	// Make sure there is no inflow
	if((2*ci-1)*res.Mass>0)
	{
		res.Mass=0;
		res.Energy=0;
		res.Momentum=(2*ci-1)*ghost.Pressure*n/abs(n);
	}
	return res;
}

Vector2D FreeFlow::CalcEdgeVelocity
(Tessellation const& /*tessellation*/,
 vector<Vector2D> const& /*point_velocities*/,
 Edge const& /*edge*/, double /*time*/) const
{
	return Vector2D(0,0);
}

bool FreeFlow::IsBoundary
(Edge const& edge,Tessellation const& tessellation)const
{
	if(edge.neighbors.first<0)
		return true;
	if(edge.neighbors.second<0)
		return true;
	return false;
}

bool FreeFlow::IsGhostCell(int i,Tessellation const& Data) const
{
	if(i>Data.GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> FreeFlow::CalcTracerFlux(Tessellation const& /*tessellation*/,
	vector<Primitive> const& /*cells*/,
	vector<vector<double> > const& tracers,double dm,
	Edge const& edge,int index,double dt,
	double /*time*/,SpatialReconstruction const& /*interp*/,
	Vector2D const& /*edge_velocity*/) const
{
  vector<double> res = tracers.at(index);
  transform(res.begin(),res.end(),res.begin(),
	    bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
  return res;
}
