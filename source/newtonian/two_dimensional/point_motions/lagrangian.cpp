#include "lagrangian.hpp"

vector<Vector2D> Lagrangian::operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
    double /*time*/, TracerStickerNames const& /*tracerstickersnames*/) const
{
  vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<res.size();++i)
    res[i] = cells[i].velocity;
  return res;
}
