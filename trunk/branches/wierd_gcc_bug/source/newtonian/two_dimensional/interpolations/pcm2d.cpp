#include "pcm2d.hpp"

void PCM2D::Prepare(Tessellation const* /*tessellation*/,
		    vector<Primitive> const& /*cells*/,
		    vector<vector<double> > const& /*tracers*/,
		    double /*dt*/,double /*time*/) {}

Primitive PCM2D::Interpolate
(Tessellation const* tessellation,
 vector<Primitive> const& cells,
 double /*dt*/,Edge const& edge, int side,
 InterpolationType interptype) const
{
  if(interptype==InBulk)
    return cells[edge.GetNeighbor(side)];
  else
    return cells[tessellation->GetOriginalIndex(edge.GetNeighbor((side+1)%2))];
}

vector<double> PCM2D::interpolateTracers
(Tessellation const* tess,vector<Primitive> const& /*cells*/,
 vector<vector<double> > const& tracers,
 double /*dt*/,
 Edge const& edge,
 int side,
 InterpolationType interp_type) const
{
  const int n = (int)tracers[edge.GetNeighbor(side)].size();
  vector<double> res(n);
  if(interp_type==InBulk)
    {
      for(int i=0;i<n;++i)
	res[i] = tracers[edge.GetNeighbor(side)][i];
    }
  else
    {
      for(int i=0;i<n;++i)
	res[i] = tracers[tess->GetOriginalIndex(edge.GetNeighbor((side+1)%2))][i];
    }
  return res;
}

