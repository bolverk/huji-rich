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

CylindricalSymmetry::CylindricalSymmetry(const Vector2D& axis):
  axis_(axis/abs(axis)) {}

double CylindricalSymmetry::calcArea(const Edge& edge_list) const
{
  const double x1 = ScalarProd(edge_list.vertices.first,axis_);
  const double x2 = ScalarProd(edge_list.vertices.second,axis_);
  const double y1 = abs(ScalarProd(edge_list.vertices.first,zcross(axis_)));
  const double y2 = abs(ScalarProd(edge_list.vertices.second,zcross(axis_)));
  return M_PI*(y1+y2)*sqrt(pow(y2-y1,2.)+pow(x2-x1,2.));
}

namespace {
  double calc_cone_segment_volume(const Vector2D& p1,
				  const Vector2D& p2,
				  const Vector2D& axis)
  {
    const double x1 = ScalarProd(p1,axis);
    const double x2 = ScalarProd(p2,axis);
    const double y1 = abs(ScalarProd(p1,zcross(axis)));
    const double y2 = abs(ScalarProd(p2,zcross(axis)));
    return (M_PI/3.)*(x2-x1)*(pow(y1,2.)+pow(y2,2.)+y1*y2);
  }

  double calc_triangular_ring_volume(const Vector2D& p1,
				     const Vector2D& p2,
				     const Vector2D& p3,
				     const Vector2D& axis)
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
