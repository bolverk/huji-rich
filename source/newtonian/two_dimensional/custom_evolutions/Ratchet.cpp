#include "Ratchet.hpp"
#include "../hydrodynamics_2d.hpp"

Ratchet::Ratchet(DIRECTION dir):
dir_(dir) {}

Conserved Ratchet::CalcFlux(Tessellation const& tess,vector<Primitive> const& cells,
	double /*dt*/,SpatialReconstruction& /*interp*/,Edge const& edge,
	Vector2D const& face_velocity,RiemannSolver const& rs,int index,
	HydroBoundaryConditions const& hbc,double /*time*/,vector<vector<double> > const& /*tracers*/)
{
	int other;
	int my_index;
	if(edge.neighbors.first==index)
	{
		my_index = 0;
		other = edge.neighbors.second;
	}
	else
	{
		my_index = 1;
		other = edge.neighbors.first;
	}
	if(hbc.IsGhostCell(other,tess))
		return Conserved();
	const Vector2D p = Parallel(edge);
	const Vector2D n = tess.GetMeshPoint(edge.neighbors.second)
		-tess.GetMeshPoint(edge.neighbors.first);
	const Vector2D outward = tess.GetMeshPoint(pair_member(edge.neighbors,1-my_index))
	  -tess.GetMeshPoint(pair_member(edge.neighbors,my_index));
	Primitive ghost = cells[other];
	if (((dir_==in)&&(ScalarProd(ghost.Velocity,outward)>0))||
		((dir_==out)&&(ScalarProd(ghost.Velocity,outward)<0)))
		ghost.Velocity = Reflect(cells[other].Velocity, p);
	Primitive left, right;
	if(0==my_index)
	{
		left = ghost;
		right = cells[other];
	}
	else
	{
		left = cells[other];
		right = ghost;
	}
	return FluxInBulk(n,p,left,right,face_velocity,rs);
}

Primitive Ratchet::UpdatePrimitive(vector<Conserved> const& /*intensives*/,
	EquationOfState const& /*eos*/,vector<Primitive>& oldcells,
	int index,Tessellation const& /*tess*/,double /*time*/,
	vector<vector<double> > const& /*tracers*/)
{
	return oldcells[index];
}


vector<double> Ratchet::UpdateTracer(int index,vector<vector<double> > const& tracers,
	vector<vector<double> > const& tracerchange,vector<Primitive> const& /*cells*/,
	Tessellation const& /*tess*/,double /*time*/)
{
	return tracers[index];
}

vector<double> Ratchet::CalcTracerFlux
	(Tessellation const& tess,vector<Primitive> const& /*cells*/,vector<vector<double> > const& tracers,
	double dm,Edge const& edge,int index,double dt,double /*time*/,
	SpatialReconstruction const& /*interp*/,Vector2D const& /*vface*/)
{
	const int other = (edge.neighbors.first==index) ?
		edge.neighbors.second : edge.neighbors.first;
	if(dir_==in)
		return vector<double>(tracers[0].size(),0);
	else
		return termwise_product(trim_top(tracers[other],tess.GetPointNo()),
		dm*dt*edge.GetLength());
}

bool Ratchet::TimeStepRelevant(void)const
{
	return false;
}
