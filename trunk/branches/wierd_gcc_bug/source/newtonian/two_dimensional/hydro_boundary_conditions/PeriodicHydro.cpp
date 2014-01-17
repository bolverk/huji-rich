#include "PeriodicHydro.hpp"

PeriodicHydro::PeriodicHydro(RiemannSolver const& rs):
rs_(rs) {}

PeriodicHydro::~PeriodicHydro(void) {}

namespace {
int calc_n(Edge const& edge,
	Tessellation const* tess)
{
	if(edge.GetNeighbor(1)>tess->GetPointNo())
		return edge.GetNeighbor(1);
	else
		return edge.GetNeighbor(0);
}
}

vector<double> PeriodicHydro::GetBoundaryTracers
	(Edge const& edge,Tessellation const* tess,
	vector<vector<double> > const& tracers,double /*time*/)const
{
	const int n = calc_n(edge,tess);
	return tracers[tess->GetOriginalIndex(n)];
}

Primitive PeriodicHydro::GetBoundaryPrimitive
	(Edge const& edge,
	Tessellation const* tess,vector<Primitive> const& cells,double /*time*/)const
{
	const int n = calc_n(edge,tess);
	return cells[tess->GetOriginalIndex(n)];
}

Conserved PeriodicHydro::CalcFlux
	(Tessellation const* tessellation,vector<Primitive> const& cells,
	Vector2D const& edge_velocity,
	Edge const& edge,SpatialReconstruction const* interp,double dt,
	double /*time*/) const
{

	const Vector2D normaldir
		(tessellation->GetMeshPoint(edge.GetNeighbor(1))-
		tessellation->GetMeshPoint(edge.GetNeighbor(0)));

	const Vector2D paraldir
		(edge.GetVertex(1) - edge.GetVertex(0));

	if(edge.GetNeighbor(0)>tessellation->GetPointNo()){
		const Primitive left =interp->Interpolate
			(tessellation,cells,dt,edge,1,Boundary);
		const Primitive right = interp->Interpolate
			(tessellation,cells,dt,edge,1,InBulk);
		return FluxInBulk
			(normaldir,paraldir,left,right,edge_velocity,rs_);
	}
	else{
		const Primitive left =interp->Interpolate
			(tessellation,cells,dt,edge,0,InBulk);
		const Primitive right = interp->Interpolate
			(tessellation,cells,dt,edge,0,Boundary);
		return FluxInBulk
			(normaldir,paraldir,left,right,edge_velocity,rs_);
	}
}

Vector2D PeriodicHydro::CalcEdgeVelocity
	(Tessellation const* tessellation,
	vector<Vector2D> const& point_velocities,
	Edge const& edge, double /*time*/) const
{
	int neigh1,neigh0;
	if(!((neigh0=edge.GetNeighbor(0))<=tessellation->GetPointNo()))
		neigh0=tessellation->GetOriginalIndex(neigh0);
	if(!((neigh1=edge.GetNeighbor(1))<=tessellation->GetPointNo()))
		neigh1=tessellation->GetOriginalIndex(neigh1);
	return tessellation->CalcFaceVelocity
		(point_velocities[neigh0],point_velocities[neigh1],
		tessellation->GetMeshPoint(edge.GetNeighbor(0)),
		tessellation->GetMeshPoint(edge.GetNeighbor(1)),
		0.5*(edge.GetVertex(0)+edge.GetVertex(1)));
}

bool PeriodicHydro::IsBoundary
	(Edge const& edge,Tessellation const* tessellation)const
{
	if((edge.GetNeighbor(0)<0)||
		(edge.GetNeighbor(0)>=tessellation->GetPointNo()))
		return true;
	if((edge.GetNeighbor(1)<0)||
		(edge.GetNeighbor(1)>=tessellation->GetPointNo()))
		return true;
	return false;
}

bool PeriodicHydro::IsGhostCell
	(int i,Tessellation const* Data) const
{
	if(i>Data->GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> PeriodicHydro::CalcTracerFlux(Tessellation const* tessellation,
	vector<Primitive> const& cells,
	vector<vector<double> > const& tracers,double dm,
	Edge const& edge,int /*index*/,double dt,
	double /*time*/,SpatialReconstruction const* interp) const
{
	vector<double> res(tracers[0].size());
	if(dm>0)
	{
		if(IsGhostCell(edge.GetNeighbor(0),tessellation))
		{
			res=interp->interpolateTracers(tessellation,cells,tracers,dt,edge,0,
				Boundary);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
		else
		{
			res=interp->interpolateTracers(tessellation,cells,tracers,dt,edge,0,
				InBulk);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
	}
	else
	{
		if(IsGhostCell(edge.GetNeighbor(0),tessellation))
		{
			res=interp->interpolateTracers(tessellation,cells,tracers,dt,edge,1,
				InBulk);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
		else
		{
			res=interp->interpolateTracers(tessellation,cells,tracers,dt,edge,1,
				Boundary);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
	}
	return res;
}
