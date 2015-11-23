#include "periodic_edge_velocities.hpp"

PeriodicEdgeVelocities::PeriodicEdgeVelocities(void) {}

vector<Vector2D> PeriodicEdgeVelocities::operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities) const
{
  const vector<Edge>& edge_list = tess.getAllEdges();
  vector<Vector2D> res(edge_list.size());
  for(size_t i=0;i<res.size();++i){
    const Edge& edge = edge_list.at(i);
    res.at(i) = tess.CalcFaceVelocity
      (point_velocities.at(static_cast<size_t>(tess.GetOriginalIndex(edge.neighbors.first))),
       point_velocities.at(static_cast<size_t>(tess.GetOriginalIndex(edge.neighbors.second))),
       tess.GetMeshPoint(edge.neighbors.first),
       tess.GetMeshPoint(edge.neighbors.second),
       calc_centroid(edge));
  }
  return res;
}
