#include "pcm2d.hpp"

bool PCM2D::WasSlopeLimited(int /*index*/)const
{
  return false;
}

void PCM2D::Prepare(Tessellation const* /*tessellation*/,
		    vector<Primitive> const& /*cells*/,
		    double /*dt*/,vector<bool> const& /*mask*/,double /*time*/) {}

Primitive PCM2D::Interpolate(Tessellation const* tessellation,
			     vector<Primitive> const& cells,
			     double /*dt*/,Edge const& edge, int side,
				 InterpolationType interptype) const
{
	if(interptype==InBulk)
		return cells[edge.GetNeighbor(side)];
	else
		return cells[tessellation->GetOriginalIndex(edge.GetNeighbor((side+1)%2))];
}

