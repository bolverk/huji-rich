#include "RigidBodyEvolve.hpp"

Conserved RigidBodyEvolve::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,	double /*dt*/,
	SpatialReconstruction& /*interpolation*/,Edge const& edge,
	Vector2D const& facevelocity,
	RiemannSolver const& rs,int index,
	HydroBoundaryConditions const& boundaryconditions,double /*time*/,
	vector<vector<double> > const& /*tracers*/)
{
	int other;
	Conserved res;
	if(edge.neighbors.first==index)
		other=edge.neighbors.second;
	else
		other=edge.neighbors.first;
	if(boundaryconditions.IsGhostCell(other,tessellation))
		return res;
	Vector2D p = Parallel(edge);
	Vector2D n = tessellation.GetMeshPoint(edge.neighbors.second)
		-tessellation.GetMeshPoint(edge.neighbors.first);
	Primitive ghost = cells[other];
	ghost.Velocity = Reflect(cells[other].Velocity, p);

	vector<Primitive> states(2);
	for(int i=0;i<2;++i)
	{
	  if(pair_member(edge.neighbors,i)==index)
			states[i] = ghost;
		else
			states[i] = cells[other];
		states[i].Velocity.Set(Projection(states[i].Velocity, n),
			Projection(states[i].Velocity, p));
	}
	res = rs.Solve(states[0], states[1],Projection(facevelocity,n));
	res.Momentum = res.Momentum.x*n/abs(n) +
		res.Momentum.y*p/abs(p);
	return res;
}

Primitive RigidBodyEvolve::UpdatePrimitive(vector<Conserved> const& /*conservedintensive*/,
	EquationOfState const& /*eos*/,vector<Primitive>& cells,int index,
	Tessellation const& /*tess*/,double /*time*/,vector<vector<double> > const& /*tracers*/)
{
	return cells[index];
}

RigidBodyEvolve::RigidBodyEvolve()
{}

RigidBodyEvolve::~RigidBodyEvolve()
{}

vector<double> RigidBodyEvolve::UpdateTracer(int index,vector<vector<double> >
	const& tracers,vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& cells,Tessellation const& tess,
	double /*time*/)
{
	return tess.GetVolume(index)*cells[index].Density*tracers[index];
}

vector<double> RigidBodyEvolve::CalcTracerFlux(Tessellation const& /*tess*/,
	vector<Primitive> const& /*cells*/,vector<vector<double> > const& tracers,
	double /*dm*/,Edge const& /*edge*/,int /*index*/,double /*dt*/,double /*time*/,
	SpatialReconstruction const& /*interp*/,Vector2D const& /*vface*/)
{
	vector<double> res(tracers[0].size(),0);
	return res;
}

bool RigidBodyEvolve::TimeStepRelevant(void)const
{
	return false;
}