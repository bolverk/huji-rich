#include "InnerRigid.hpp"


InnerRigid::InnerRigid(int n,RiemannSolver const& rs,HydroBoundaryConditions *hbc): 
PointNum(n),rs_(rs),hbc_(hbc),rhbc(rs)
{}

InnerRigid::~InnerRigid(){}

bool InnerRigid::IsGhostCell(int i,Tessellation const* tess)const
{
	if(i<PointNum)
		return true;
	else
		return hbc_->IsGhostCell(i,tess);
}

bool InnerRigid::IsBoundary(Edge const& edge,Tessellation const* tess)const
{
	if(edge.GetNeighbor(0)<PointNum)
		return true;
	if(edge.GetNeighbor(1)<PointNum)
		return true;
	return hbc_->IsBoundary(edge,tess);
}

Vector2D InnerRigid::CalcEdgeVelocity(Tessellation const* tessellation,
	vector<Vector2D> const& point_velocities,Edge const& edge,double time) const
{
	if(hbc_->IsBoundary(edge,tessellation))
		return hbc_->CalcEdgeVelocity(tessellation,point_velocities,edge,
		time);
	int neigh1=edge.GetNeighbor(1),neigh0=edge.GetNeighbor(0);	
	return tessellation->CalcFaceVelocity(point_velocities[neigh0],point_velocities[neigh1],
		tessellation->GetMeshPoint(edge.GetNeighbor(0)),
		tessellation->GetMeshPoint(edge.GetNeighbor(1)),
		0.5*(edge.GetVertex(0)+edge.GetVertex(1)));
}

Conserved InnerRigid::CalcFlux(Tessellation const* tessellation,
	vector<Primitive> const& cells,Vector2D const& edge_velocity,
	Edge const& edge,SpatialReconstruction const* interp,double dt,
	double time) const
{
	if(hbc_->IsBoundary(edge,tessellation))
		return hbc_->CalcFlux(tessellation,cells,edge_velocity,edge,interp,
		dt,time);
	else
	{
		int ci;
		if((edge.GetNeighbor(0)<PointNum)&&(edge.GetNeighbor(1)>=PointNum))
			ci = 1;
		else
		{
			if((edge.GetNeighbor(1)<PointNum)&&(edge.GetNeighbor(0)>=
				PointNum))
				ci = 0;
			else
			{
				if((edge.GetNeighbor(1)<PointNum)&&(edge.GetNeighbor(0)<
					PointNum))
					return Conserved();
				else
					throw UniversalError("Error in InnerRigid");
			}
		}
		return CalcFluxCi(tessellation,cells,edge_velocity,&rs_,edge,
			interp,dt,ci);
	}
}

Primitive InnerRigid::GetBoundaryPrimitive(Edge const& edge,
	  Tessellation const* tess,vector<Primitive> const& cells,
	  double time)const
{
	if(hbc_->IsBoundary(edge,tess))
		return hbc_->GetBoundaryPrimitive(edge,tess,cells,time);
	else
	{
		const Vector2D p = Parallel(edge);
		int ci;
		if(edge.GetNeighbor(0)<PointNum)
			ci=0;
		else
			ci=1;
		Primitive res =cells[edge.GetNeighbor(ci)];
		res.Velocity = Reflect(res.Velocity,p);
		return res;
	}
}

vector<double> InnerRigid::GetBoundaryTracers(Edge const& edge,
	  Tessellation const* tess,vector<vector<double> > const& tracers,
	  double time)const
{
	if(hbc_->IsBoundary(edge,tess))
		return hbc_->GetBoundaryTracers(edge,tess,tracers,time);
	else
	{
		if(edge.GetNeighbor(0)<PointNum)
			return tracers[edge.GetNeighbor(0)];
		else
			return tracers[edge.GetNeighbor(1)];
	}
}

vector<double> InnerRigid::CalcTracerFlux(Tessellation const* tessellation,
	  vector<vector<double> > const& tracers,double dm,
	  Edge const& edge,int index,double dt,
	  double time,SpatialReconstruction const* interp)const
{
	if(hbc_->IsBoundary(edge,tessellation))
		return hbc_->CalcTracerFlux(tessellation,tracers,dm,edge,index,
		dt,time,interp);
	else
	{
		vector<double> res(tracers[0].size(),0);
		return res;
	}
}

Conserved InnerRigid::CalcFluxCi(Tessellation const* tessellation,
	vector<Primitive> const& cells,Vector2D const& edge_velocity,
	RiemannSolver const* rs,Edge const& edge,
	SpatialReconstruction const* interp,double dt,int ci) const
{
	const Vector2D p = Parallel(edge);
	const Vector2D n = Normal(edge, tessellation);
	Primitive ghost = interp->Interpolate
		(tessellation,cells,dt,edge,ci,
		InBulk);
	const Primitive othercell=ghost;
	ghost.Velocity = 2*edge_velocity+Reflect(ghost.Velocity, p);
	vector<Primitive> states(2);
	for(int i=0;i<2;i++)
	{
		if(IsGhostCell(edge.GetNeighbor(i),tessellation))
			states[i] = ghost;
		else
			states[i] = othercell;
		states[i].Velocity.Set
			(Projection(states[i].Velocity, n),
			Projection(states[i].Velocity, p));
	}
	Conserved res = rs->Solve(states[0], states[1],Projection(edge_velocity,n));
	res.Momentum = res.Momentum.x*n/abs(n) +
		res.Momentum.y*p/abs(p);
	return res;
}

int InnerRigid::GetPointNum(void)const
{
	return PointNum;
}
