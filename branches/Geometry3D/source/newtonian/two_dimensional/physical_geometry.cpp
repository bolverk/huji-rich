#include "physical_geometry.hpp"
#include "../../tessellation/triangle_area.hpp"

PhysicalGeometry::~PhysicalGeometry(void) {}

SlabSymmetry::SlabSymmetry(void) {}

double SlabSymmetry::calcArea(const Edge& edge) const
{
  return abs(edge.vertices.first-edge.vertices.second);
}

double SlabSymmetry::calcVolume(const vector<Edge>& edge_list) const
{
  double res = 0;
  const Vector2D anchor = edge_list[0].vertices.first;
  for(size_t i=1;i<edge_list.size();++i)
    res += calc_triangle_area(anchor,edge_list[i].vertices.first,
			      edge_list[i].vertices.second);
  return res;
}

Axis::Axis(const Vector2D& origin_i,
	   const Vector2D& direction_i):
  origin(origin_i), direction(direction_i/abs(direction_i)) {}

namespace {
  Vector2D change_coordinate(const Vector2D& v,
			     const Axis& axis)
  {
    return Vector2D(ScalarProd(v-axis.origin,
			       axis.direction),
		    std::abs(ScalarProd(v-axis.origin,
					zcross(axis.direction))));
  }
}

CylindricalSymmetry::CylindricalSymmetry
(const Vector2D& origin, const Vector2D& direction):
  axis_(origin,direction) {}

double CylindricalSymmetry::calcArea(const Edge& edge) const
{
  const Vector2D p1 = change_coordinate(edge.vertices.first, axis_);
  const Vector2D p2 = change_coordinate(edge.vertices.second, axis_);
  return M_PI*(p1.y+p2.y)*sqrt(pow(p2.y-p1.y,2.)+pow(p2.x-p1.x,2.));
}

namespace {
  double calc_cone_segment_volume(const Vector2D& p1,
				  const Vector2D& p2,
				  const Axis& axis)
  {
    const Vector2D q1 = change_coordinate(p1, axis);
    const Vector2D q2 = change_coordinate(p2, axis);
    return (M_PI/3.)*(q2.x-q1.x)*(pow(q1.y,2.)+pow(q2.y,2.)+q1.y*q2.y);
  }

  double calc_triangular_ring_volume(const Vector2D& p1,
				     const Vector2D& p2,
				     const Vector2D& p3,
				     const Axis& axis)
  {
    return std::abs(calc_cone_segment_volume(p1,p2,axis)+
		    calc_cone_segment_volume(p2,p3,axis)+
		    calc_cone_segment_volume(p3,p1,axis));
  }
}

double CylindricalSymmetry::calcVolume(const vector<Edge>& edge_list) const
{
  double res = 0;
  const Vector2D anchor = edge_list[0].vertices.first;
  for(size_t i=1;i<edge_list.size();++i)
    res += calc_triangular_ring_volume(anchor,
				       edge_list[i].vertices.first,
				       edge_list[i].vertices.second,
				       axis_);
  return res;
}

const Axis& CylindricalSymmetry::getAxis(void) const
{
  return axis_;
}
