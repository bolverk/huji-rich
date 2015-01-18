#include "RigidWallHydro.hpp"

RigidWallHydro::RigidWallHydro(RiemannSolver const& rs):
rs_(rs) {}

RigidWallHydro::~RigidWallHydro(void) {}

vector<double> RigidWallHydro::GetBoundaryTracers(Edge const& edge,Tessellation const& /*Data*/,
	vector<vector<double> > const& tracers,double /*time*/)const
{
	if(edge.neighbors.first==-1)
	  return tracers[static_cast<size_t>(edge.neighbors.second)];
	else
	  return tracers[static_cast<size_t>(edge.neighbors.first)];
}

Primitive RigidWallHydro::GetBoundaryPrimitive(Edge const& edge,
	Tessellation const& /*Data*/,vector<Primitive> const& cells,
	double /*time*/)const
{
	const Vector2D p = Parallel(edge);
	Primitive res =cells[static_cast<size_t>(pair_member(edge.neighbors,GetRealCell(edge)))];
	res.Velocity = Reflect(res.Velocity,p);
	return res;
}

int RigidWallHydro::GetRealCell(Edge const& edge)const
{
	if((edge.neighbors.first==-1)&&(edge.neighbors.second!=-1))
		return 1;
	else{
		if((edge.neighbors.second==-1)&&(edge.neighbors.first!=-1))
			return 0;
		else{
		  throw UniversalError("Boundary condition called for bulk cell");
		}
	}
}

Conserved RigidWallHydro::CalcFluxCi
	(Tessellation const& tessellation,
	vector<Primitive> const& cells,Vector2D const& edge_velocity,
	RiemannSolver const& rs,Edge const& edge,
	SpatialReconstruction const& interp,double dt,int ci) const
{
	const Vector2D p = Parallel(edge);
	const Vector2D n = Normal(edge, tessellation);
	Primitive ghost = interp.Interpolate
		(tessellation,cells,dt,edge,ci,
		InBulk,edge_velocity);
	const Primitive othercell=ghost;
	ghost.Velocity = Reflect(ghost.Velocity, p);
	vector<Primitive> states(2);
	for(int i=0;i<2;i++){
	  if(IsGhostCell(pair_member(edge.neighbors,i),tessellation))
	    states[static_cast<size_t>(i)] = ghost;
		else
		  states[static_cast<size_t>(i)] = othercell;
	  states[static_cast<size_t>(i)].Velocity.Set
	    (Projection(states[static_cast<size_t>(i)].Velocity, n),
	     Projection(states[static_cast<size_t>(i)].Velocity, p));
	}
	Conserved res = rs.Solve(states[0], states[1],Projection(edge_velocity,n));
	res.Momentum = res.Momentum.x*n/abs(n) +
		res.Momentum.y*p/abs(p);
	return res;
}

Conserved RigidWallHydro::CalcFlux
	(Tessellation const& tessellation,
	vector<Primitive> const& cells,Vector2D const& edge_velocity,
	Edge const& edge,
	SpatialReconstruction const& interp,
	double dt, double /*time*/) const
{
	const int ci=GetRealCell(edge);
	return CalcFluxCi(tessellation,cells,edge_velocity,rs_,edge,interp,dt,ci);
}

Vector2D RigidWallHydro::CalcEdgeVelocity
	(Tessellation const& /*tessellation*/,
	vector<Vector2D> const& /*point_velocities*/,
	Edge const& /*edge*/,
	double /*time*/) const
{
	return Vector2D(0,0);
}

bool RigidWallHydro::IsBoundary(Edge const& edge,Tessellation const& /*Data*/)const
{
	if(edge.neighbors.first<0)
		return true;
	if(edge.neighbors.second<0)
		return true;
	return false;
}

bool RigidWallHydro::IsGhostCell(int i,Tessellation const& Data)const
{
	if(i>Data.GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> RigidWallHydro::CalcTracerFlux(Tessellation const& /*tessellation*/,
	vector<Primitive> const& /*cells*/,
		vector<vector<double> > const& tracers,double /*dm*/,
		Edge const& /*edge*/,int /*index*/,double /*dt*/,
		double /*time*/,SpatialReconstruction const& /*interp*/,
		Vector2D const& /*edge_velocity*/) const
{
  vector<double> res(tracers[0].size(),0);
  return res;
}
