#include "Ratchet.hpp"
#include "../hydrodynamics_2d.hpp"

Ratchet::Ratchet(DIRECTION dir, const PhysicalGeometry& pg):
  dir_(dir),cell_(Primitive(1,1,Vector2D(0,0),1,1)), pg_(pg) {}

Conserved Ratchet::CalcFlux(Tessellation const& tess,
			    vector<Primitive> const& cells,
			    double /*dt*/,
			    SpatialReconstruction& /*interp*/,
			    Edge const& edge,
			    Vector2D const& face_velocity,
			    RiemannSolver const& rs,int index,
			    HydroBoundaryConditions const& hbc,
			    double /*time*/,
			    vector<vector<double> > const& /*tracers*/)
{
  const int other = (edge.neighbors.first==index) ?
    edge.neighbors.second : edge.neighbors.first;
  if(hbc.IsGhostCell(other,tess))
    return Conserved();
  const int my_index = (edge.neighbors.first==index) ? 0 : 1;
  const Vector2D p = Parallel(edge);
  const Vector2D n = tess.GetMeshPoint(edge.neighbors.second)
    -tess.GetMeshPoint(edge.neighbors.first);
  const Vector2D outward =
    tess.GetMeshPoint(pair_member(edge.neighbors,1-my_index))
    - tess.GetMeshPoint(pair_member(edge.neighbors,my_index));
  Primitive ghost = cells[static_cast<size_t>(other)];
  if (((dir_==in)&&(ScalarProd(ghost.Velocity,outward)>0))||
      ((dir_==out)&&(ScalarProd(ghost.Velocity,outward)<0)))
    ghost.Velocity = Reflect(cells[static_cast<size_t>(other)].Velocity, p);
  Primitive left, right;
  if(0==my_index)
    {
      left = ghost;
      right = cells[static_cast<size_t>(other)];
    }
  else
    {
      left = cells[static_cast<size_t>(other)];
      right = ghost;
    }
  return FluxInBulk(n,p,left,right,face_velocity,rs);
}

Primitive Ratchet::UpdatePrimitive(vector<Conserved> const& /*intensives*/,
				   EquationOfState const& /*eos*/,vector<Primitive>& /*oldcells*/,
				   int /*index*/,Tessellation const& /*tess*/,double /*time*/,
	vector<vector<double> > const& /*tracers*/)
{
	return cell_;
}

vector<double> Ratchet::UpdateTracer(int index,vector<vector<double> > const& tracers,
				     vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& /*cells*/,
	Tessellation const& /*tess*/,double /*time*/)
{
  return tracers[static_cast<size_t>(index)];
}

vector<double> Ratchet::CalcTracerFlux
(Tessellation const& /*tess*/,vector<Primitive> const& /*cells*/,vector<vector<double> > const& tracers,
	double dm,Edge const& edge,int index,double dt,double /*time*/,
	SpatialReconstruction const& /*interp*/,Vector2D const& /*vface*/)
{
	const int other = (edge.neighbors.first==index) ?
	  edge.neighbors.second : edge.neighbors.first;
	return termwise_product(tracers[static_cast<size_t>(other)],
				dm*dt*pg_.calcArea(edge));
}

bool Ratchet::TimeStepRelevant(void)const
{
	return false;
}

bool Ratchet::isRelevantToInterpolation(void)const
{
	return false;
}
