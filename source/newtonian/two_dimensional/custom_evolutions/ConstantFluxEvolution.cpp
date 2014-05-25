#include "ConstantFluxEvolution.hpp"

ConstantFluxEvolution::ConstantFluxEvolution(Primitive const& cell,
	vector<double> const& tracer,EquationOfState const& eos,bool entropycalc)
:cell_(cell),tracer_(tracer),eos_(eos),entropy_(entropycalc){}

ConstantFluxEvolution::~ConstantFluxEvolution(void)
{}

Conserved ConstantFluxEvolution::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,	double dt,
	SpatialReconstruction& interpolation,Edge const& edge,
	Vector2D const& facevelocity,RiemannSolver const& /*rs*/,int index,
	HydroBoundaryConditions const& bc,double time,
	vector<vector<double> > const& /*tracers*/)
{
	if(bc.IsBoundary(edge,tessellation))
		return bc.CalcFlux(tessellation,cells,facevelocity,edge,interpolation,dt,
		time);
	else
	{
		Vector2D normaldir = tessellation.GetMeshPoint(edge.GetNeighbor(1))-
			tessellation.GetMeshPoint(edge.GetNeighbor(0));
		normaldir=normaldir/abs(normaldir);
		double v=abs(cell_.Velocity);
		Primitive cell(cell_);
		if(edge.GetNeighbor(1)==index)
			cell.Velocity=-v*normaldir;
		else
			cell.Velocity=v*normaldir;
		return Primitive2Flux(cell,normaldir);
	}
}

Primitive ConstantFluxEvolution::UpdatePrimitive
	(vector<Conserved> const& /*conservedintensive*/,
	EquationOfState const& /*eos*/,
	vector<Primitive>& cells,int index,Tessellation const& /*tess*/,
	double /*time*/,vector<vector<double> > const& /*tracers*/)
{
	Primitive res=cells[index];
	return res;
}

vector<double> ConstantFluxEvolution::UpdateTracer(int index,vector<vector<double> >
	const& /*tracers*/,vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& cells,
	Tessellation const& tess,double /*time*/)
{
	if(entropy_)
	{
		vector<double> res(tracer_);
		res[0]=eos_.dp2s(cells[index].Density,cells[index].Pressure);
		return tess.GetVolume(index)*cells[index].Density*res;
	}
	else
		return tess.GetVolume(index)*cells[index].Density*tracer_;
}

vector<double> ConstantFluxEvolution::CalcTracerFlux(Tessellation const& /*tess*/,
	vector<Primitive> const& /*cells*/,vector<vector<double> > const& /*tracers*/,
	double dm,Edge const& edge,int /*index*/,double dt,double /*time*/,
	SpatialReconstruction const& /*interp*/,Vector2D const& /*vface*/)
{
	vector<double> res(tracer_);
	if(entropy_)
		res[0]=eos_.dp2s(cell_.Density,cell_.Pressure);
	transform(res.begin(),res.end(),res.begin(),
		bind1st(multiplies<double>(),abs(dm)*dt*edge.GetLength()));	
	return res;
}
