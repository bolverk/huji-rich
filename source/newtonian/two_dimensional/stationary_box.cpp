#include "stationary_box.hpp"

StationaryBox::StationaryBox(void) {}

namespace {
  bool is_outer_edge
  (const Tessellation& tess,
   const Edge& edge)
  {
    const int n1 =
      tess.GetOriginalIndex
       (edge.neighbors.first);
    const int n2 =
      tess.GetOriginalIndex
       (edge.neighbors.second);
    return n1==n2;
  }
}

vector<Vector2D> StationaryBox::operator()
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities) const
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  vector<Vector2D> res(edge_list.size());
  const size_t nloop = res.size();
  for(size_t i=0;i<nloop;++i)
  {
    const Edge& edge = edge_list[i];
    res[i] = is_outer_edge(tess,edge) ?
      Vector2D() :
      tess.CalcFaceVelocity
      (point_velocities.at(static_cast<size_t>(edge.neighbors.first)),
       point_velocities.at(static_cast<size_t>(edge.neighbors.second)),
       tess.GetMeshPoint(edge.neighbors.first),
       tess.GetMeshPoint(edge.neighbors.second),
       calc_centroid(edge));
  }
  return res;
}
