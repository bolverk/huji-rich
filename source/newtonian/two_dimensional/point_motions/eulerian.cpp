#include "eulerian.hpp"

vector<Vector2D> Eulerian::operator()
(Tessellation const& tessellation,
 vector<ComputationalCell> const& /*primitives*/,double /*time*/, TracerStickerNames const& /*tracerstickersnames*/) const
{
  return vector<Vector2D>
    (static_cast<size_t>(tessellation.GetPointNo()));
}
