#include "noh_hbc.hpp"

using std::multiplies;

NohHBC::NohHBC(Vector2D const& center,double d0,double v0,double p0 ):
  center_(center), d0_(d0), v0_(v0),p0_(p0) {}

Conserved NohHBC::CalcFlux(Tessellation const& tess,vector<Primitive> const& /*cells*/,
	Vector2D const& /*edge_velocity*/,Edge const& edge,
	SpatialReconstruction const& /*interp*/,double /*dt*/,double time) const
{
	const Vector2D edge_cen = 0.5*(edge.vertices.first+
				       edge.vertices.second);
	const double r = abs(center_-edge_cen);
	if(r==0)
		throw UniversalError("Error in NohHBC. Radius is zero");
	const double d = d0_*(1+time*v0_/r);
	Primitive p;
	p.Density = d;
	p.Velocity = v0_*(center_-edge_cen)/r;
	return Primitive2Flux(p,Normal(edge,tess));
}

Primitive NohHBC::GetBoundaryPrimitive(Edge const& edge,
	Tessellation const& /*tess*/,vector<Primitive> const& /*cells*/,double
	time)const
{
	const Vector2D edge_cen = 0.5*(edge.vertices.first+
				       edge.vertices.second);
	const double r = abs(center_-edge_cen);
	if(r==0)
		throw UniversalError("Error in NohHBC. Radius is zero");
	const double d = d0_*(1+time*v0_/r);
	Primitive p;
	p.Density = d;
	p.Velocity = v0_*(center_-edge_cen)/r;
	p.Pressure=p0_;
	return p;
}

Vector2D NohHBC::CalcEdgeVelocity
	(Tessellation const& /*tessellation*/,
	vector<Vector2D> const& /*point_velocities*/,
	Edge const& /*edge*/, double /*time*/) const
{
	Vector2D res(0,0);
	return res;
}

bool NohHBC::IsBoundary(Edge const& edge,Tessellation const& /*tessellation*/)const
{
	if(edge.neighbors.first<0)
		return true;
	if(edge.neighbors.second<0)
		return true;
	return false;
}

bool NohHBC::IsGhostCell(int i,Tessellation const& tess) const
{
	if(i>tess.GetPointNo()||i<0)
		return true;
	else
		return false;
}

vector<double> NohHBC::GetBoundaryTracers(Edge const& edge,
	Tessellation const& tess,vector<vector<double> > const& tracers,double
	time)const
{
	vector<double> res;
	if(tracers.empty())
		throw UniversalError("Noh problem must use cold flows flag");
	const Vector2D edge_cen = 0.5*(edge.vertices.first+
				       edge.vertices.second);
	const double r = abs(center_-edge_cen);
	if(r==0)
		throw UniversalError("Error in NohHBC. Radius is zero");
	const double d = d0_*(1+time*v0_/r);
	Primitive p;
	p.Density = d;
	p.Pressure=p0_;
	res.push_back(p.Pressure*pow(d,-1.66666));
	for(size_t i=1;i<tracers[0].size();++i)
	{
		if(IsGhostCell(edge.neighbors.second,tess))
			res.push_back(tracers[edge.neighbors.first][i]);
		else
			res.push_back(tracers[edge.neighbors.second][i]);
	}

	return res;
}

vector<double> NohHBC::CalcTracerFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells,
	  vector<vector<double> > const& tracers,double dm,
	  Edge const& edge,int /*index*/,double dt,
	  double time,SpatialReconstruction const& interp,
	  Vector2D const& vface) const
{
	vector<double> res;
	if(IsGhostCell(edge.neighbors.first,tessellation))
	{
		if(dm>0)
		{
			res=GetBoundaryTracers(edge,tessellation,tracers,time);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
		else
		{
			res=interp.interpolateTracers(tessellation,cells,tracers,dt,edge,1,
				InBulk,vface);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
	}
	else
	{
		if(dm<0)
		{
			res=GetBoundaryTracers(edge,tessellation,tracers,time);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
		else
		{
			res=interp.interpolateTracers(tessellation,cells,tracers,dt,edge,0,
				InBulk,vface);
			transform(res.begin(),res.end(),res.begin(),
				bind1st(multiplies<double>(),dm*dt*edge.GetLength()));
		}
	}
	return res;
}
