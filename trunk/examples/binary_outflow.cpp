#include "binary_outflow.hpp"

BinaryOutflow::BinaryOutflow(double mdot,double e,double cs,int nx,int ny):
mdot_(mdot),e_(e),cs_(cs),nx_(nx),ny_(ny)
{}

int BinaryOutflow::GetCellIndex(void)
{
	return cell_index_;
}

void BinaryOutflow::Init(Tessellation const& tess,double xc,double yc)
{	
	// Find the outflow cell index
	cell_index_=0;
	double distance=tess.GetMeshPoint(0).distance(Vector2D(xc,yc));
	for(int i=1;i<nx_*ny_;++i)
	{
		if(tess.GetMeshPoint(i).distance(Vector2D(xc,yc))<distance)
		{
			distance=tess.GetMeshPoint(i).distance(Vector2D(xc,yc));
			cell_index_=i;
		}
	}
}

Conserved BinaryOutflow::Calculate(Tessellation const* tess,
	vector<Primitive> const& /*cells*/,int point,vector<Conserved> const& /*fluxes*/,
	vector<Vector2D> const& /*point_velocity*/,HydroBoundaryConditions const* /*hbc*/,
	vector<vector<double> > const& /*tracers*/,vector<double> & /*dtracer*/,
	double t,double dt)
{
	Conserved res;
	// Are we close to the source?? We want to add to the 8 closest cells
	if((point==(cell_index_+1))||(point==(cell_index_-1))||(point==(cell_index_+nx_))
		||(point==(cell_index_-nx_))||(point==(cell_index_+1+nx_))||
		(point==(cell_index_+1-nx_))||(point==(cell_index_+nx_-1))||
		(point==(cell_index_-1-nx_)))
	{
		res.Mass=mdot_;
		res.Energy=mdot_*e_;
		Vector2D direction=tess->GetMeshPoint(point)-tess->GetMeshPoint(cell_index_);
		direction=direction/abs(direction);
		res.Momentum=mdot_*cs_*direction;
	}
	return res;
}