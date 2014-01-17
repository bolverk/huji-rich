#include "pcm_scalar.hpp"

void PCMScalar::Prepare(Tessellation const* /*tess*/,
			vector<vector<double> > const& /*tracers*/,
			double /*dt*/,
			double /*time*/) {}

vector<double> PCMScalar::Interpolate(Tessellation const* tess,
					vector<vector<double> > const& tracers,
					double /*dt*/,
					Edge const& edge,
					int side,
					SInterpolationType interptype) const
{
  if(interptype==SInBulk)
    return tracers[edge.GetNeighbor(side)];
  else
    return tracers[tess->GetOriginalIndex(edge.GetNeighbor((side+1)%2))];
}
