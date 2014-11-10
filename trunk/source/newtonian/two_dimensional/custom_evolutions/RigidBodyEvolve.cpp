#include "RigidBodyEvolve.hpp"
#include "../hydrodynamics_2d.hpp"

namespace {

  Primitive reflect(const Primitive& cell,
		    const Vector2D& axis)
  {
    Primitive res = cell;
    res.Velocity = Reflect(cell.Velocity, axis);
    return res;
  }

  Primitive rotate(const Primitive& cell,
		   const Vector2D& n,
		   const Vector2D& p)
  {
    Primitive res = cell;
    res.Velocity = Vector2D(Projection(res.Velocity, n),
			    Projection(res.Velocity, p));
    return res;
  }

  Conserved rotate(const Conserved& c,
		   const Vector2D& n,
		   const Vector2D& p)
  {
    Conserved res = c;
    res.Momentum =  res.Momentum.x*n/abs(n) +
      res.Momentum.y*p/abs(p);
    return res;
  }
}

Conserved RigidBodyEvolve::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,	double /*dt*/,
	SpatialReconstruction& /*interpolation*/,Edge const& edge,
	Vector2D const& facevelocity,
	RiemannSolver const& rs,int index,
	HydroBoundaryConditions const& boundaryconditions,double /*time*/,
	vector<vector<double> > const& /*tracers*/)
{
  const int other = get_other_index(edge,index);
  if(boundaryconditions.IsGhostCell(other,tessellation))
    return Conserved();
  const Vector2D p = Parallel(edge);
  const Vector2D n = tessellation.GetMeshPoint(edge.neighbors.second)
    -tessellation.GetMeshPoint(edge.neighbors.first);
  const Primitive ghost = reflect(cells[(size_t)other],p);
  const std::pair<Primitive, Primitive> states(rotate(edge.neighbors.first==index ?
						      ghost : cells[(size_t)other],n,p),
					       rotate(edge.neighbors.second==index ?
						      ghost : cells[(size_t)other],n,p));
  return rotate(rs.Solve(states.first, states.second, Projection(facevelocity,n)),n,p);
}

Primitive RigidBodyEvolve::UpdatePrimitive(vector<Conserved> const& /*conservedintensive*/,
					   EquationOfState const& /*eos*/,vector<Primitive>& /*cells*/,int /*index*/,
	Tessellation const& /*tess*/,double /*time*/,vector<vector<double> > const& /*tracers*/)
{
	return Primitive(1,1,Vector2D(0,0),1,1);
}

RigidBodyEvolve::RigidBodyEvolve()
{}

RigidBodyEvolve::~RigidBodyEvolve()
{}

vector<double> RigidBodyEvolve::UpdateTracer(int index,vector<vector<double> >
	const& tracers,vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& cells,Tessellation const& tess,
	double /*time*/)
{
  return tess.GetVolume(index)*cells[(size_t)index].Density*tracers[(size_t)index];
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
