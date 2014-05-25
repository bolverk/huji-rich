#include "RigidBodyEvolve.hpp"

Conserved RigidBodyEvolve::CalcFlux(Tessellation const* tessellation,
	vector<Primitive> const& cells,	double /*dt*/,
	SpatialReconstruction* /*interpolation*/,Edge const& edge,
	Vector2D const& facevelocity,
	RiemannSolver const& rs,int index,
	HydroBoundaryConditions const* boundaryconditions,double /*time*/,
	vector<vector<double> > const& /*tracers*/)
{
	int other;
	Conserved res;
	if(edge.GetNeighbor(0)==index)
		other=edge.GetNeighbor(1);
	else
		other=edge.GetNeighbor(0);
	if(boundaryconditions->IsGhostCell(other,tessellation))
		return res;
	Vector2D p = Parallel(edge);
	Vector2D n = tessellation->GetMeshPoint(edge.GetNeighbor(1))
		-tessellation->GetMeshPoint(edge.GetNeighbor(0));
	Primitive ghost = cells[other];
	ghost.Velocity = Reflect(cells[other].Velocity, p);

	vector<Primitive> states(2);
	for(int i=0;i<2;++i)
	{
		if(edge.GetNeighbor(i)==index)
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
	EquationOfState const* /*eos*/,vector<Primitive>& cells,int index)
{
	return cells[index];
}


RigidBodyEvolve::RigidBodyEvolve()
{}
