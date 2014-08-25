#include "pcm2d.hpp"
#include "../../../misc/utils.hpp"

PCM2D::PCM2D(void): slopes_(vector<ReducedPrimitiveGradient2D> ()){}

void PCM2D::Prepare(Tessellation const& /*tessellation*/,
		    vector<Primitive> const& /*cells*/,
		    vector<vector<double> > const& /*tracers*/,
		    double /*dt*/,double /*time*/) {}

Primitive PCM2D::Interpolate
(Tessellation const& tessellation,
 vector<Primitive> const& cells,
 double /*dt*/,Edge const& edge, int side,
 InterpolationType interptype,Vector2D const& /*vface*/) const
{
  if(interptype==InBulk)
    return cells[pair_member(edge.neighbors,side)];
  else
    return cells[pair_member(edge.neighbors,(side+1)%2)];
}

vector<ReducedPrimitiveGradient2D>& PCM2D::GetGradients(void)
{
	return slopes_;
}

vector<double> PCM2D::interpolateTracers
(Tessellation const& tess,vector<Primitive> const& /*cells*/,
 vector<vector<double> > const& tracers,
 double /*dt*/,
 Edge const& edge,
 int side,
 InterpolationType interp_type,Vector2D const& /*vface*/) const
{
  const int n = (int)tracers[pair_member(edge.neighbors,side)].size();
  vector<double> res(n);
  if(interp_type==InBulk)
    {
      for(int i=0;i<n;++i)
	res[i] = tracers[pair_member(edge.neighbors,side)][i];
    }
  else
    {
      for(int i=0;i<n;++i)
	res[i] = tracers[pair_member(edge.neighbors,(side+1)%2)][i];
    }
  return res;
}