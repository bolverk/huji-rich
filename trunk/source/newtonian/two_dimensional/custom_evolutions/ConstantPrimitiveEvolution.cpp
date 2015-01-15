#include "ConstantPrimitiveEvolution.hpp"

using std::multiplies;

bool ConstantPrimitiveEvolution::isRelevantToInterpolation(void) const
{
	return true;
}

bool ConstantPrimitiveEvolution::ShouldForceTracerReset(void)const
{
	return true;
}

ConstantPrimitiveEvolution::ConstantPrimitiveEvolution(Primitive const& prim,
	vector<double> const& tracer,bool mass_count,int n): prim_(prim),tracer_(tracer),
	mass_count_(mass_count),mass_flux(0),mass_fluxt(0),N_(n)
{}

ConstantPrimitiveEvolution::~ConstantPrimitiveEvolution(void)
{}

Conserved ConstantPrimitiveEvolution::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,	double dt,
	SpatialReconstruction& interpolation,Edge const& edge,
	Vector2D const& facevelocity,RiemannSolver const& rs,int index,
	HydroBoundaryConditions const& bc,double time,vector<vector<double> > const& tracers)
{
	if(bc.IsBoundary(edge,tessellation))
		return bc.CalcFlux(tessellation,cells,facevelocity,edge,interpolation,dt,time);
	else
	{
		Vector2D normaldir = tessellation.GetMeshPoint(edge.neighbors.second)-
			tessellation.GetMeshPoint(edge.neighbors.first);

		Vector2D paraldir = edge.vertices.second - edge.vertices.first;

		Primitive left = interpolation.Interpolate
			(tessellation, cells, dt, edge, 0,InBulk,facevelocity);
		Primitive right = interpolation.Interpolate
			(tessellation, cells, dt, edge, 1,InBulk,facevelocity);
		Conserved res(FluxInBulk(normaldir,paraldir,left,right,facevelocity,rs));
		int n0=edge.neighbors.first;
		int n1=edge.neighbors.second;
		// Do not allow outflow from this region
		if(((n0<N_&&res.Mass>0)||(n1<N_&&res.Mass<0))&&(n0>=N_||n1>=N_))
		{
		return Conserved();
		}
		if(mass_count_)
		{
		  if((n0>=N_&&tracers[static_cast<size_t>(n0)][1]>1e-6)||(n1>=N_&&tracers[static_cast<size_t>(n1)][1]>1e-6))
			{
				if(n0==index)
				{
					mass_flux+=res.Mass*edge.GetLength()*dt;
					mass_fluxt+=res.Mass*edge.GetLength()*dt*tracers[static_cast<size_t>(n1)][1];
				}
				else
				{
					mass_flux-=res.Mass*edge.GetLength()*dt;
					mass_fluxt-=res.Mass*edge.GetLength()*dt*tracers[static_cast<size_t>(n0)][1];
				}
			}
		}
		return res;
	}
}

Primitive ConstantPrimitiveEvolution::UpdatePrimitive
	(vector<Conserved> const& /*conservedintensive*/,
	EquationOfState const& /*eos*/,
	vector<Primitive>& /*cells*/,int /*index*/,Tessellation const& /*tess*/,
	double /*time*/,vector<vector<double> > const& /*tracers*/)
{
	return prim_;
}

void ConstantPrimitiveEvolution::SetMassFlux(double m)
{
	mass_flux=m;
}

double ConstantPrimitiveEvolution::GetMassFlux(void) const
{
	return mass_flux;
}

vector<double> ConstantPrimitiveEvolution::UpdateTracer
(int index,vector<vector<double> >
 const& /*tracers*/,vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& /*cells*/,
 Tessellation const& tess,double /*time*/)
{
  return prim_.Density*tess.GetVolume(index)*tracer_;
}

vector<double> ConstantPrimitiveEvolution::CalcTracerFlux
(Tessellation const& tess,
 vector<Primitive> const& cells,vector<vector<double> > const& tracers,
 double dm,Edge const& edge,int /*index*/,double dt,double /*time*/,
 SpatialReconstruction const& interp,Vector2D const& vface)
{
	vector<double> res(tracers[0].size());
	if(dm>0)
	{
		res=interp.interpolateTracers(tess,cells,tracers,dt,edge,0,
			InBulk,vface);
		transform(res.begin(),res.end(),res.begin(),
			bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
	}
	else
	{
		res=interp.interpolateTracers(tess,cells,tracers,dt,edge,1,
			Boundary,vface);
		transform(res.begin(),res.end(),res.begin(),
			bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
	}
	return res;
}

bool ConstantPrimitiveEvolution::TimeStepRelevant(void)const
{
	return true;
}
