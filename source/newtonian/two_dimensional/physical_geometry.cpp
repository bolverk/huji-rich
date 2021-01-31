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

Vector2D SlabSymmetry::calcCentroid(const vector<Vector2D>& chull) const
{
	Vector2D res;
	double V = 0;
	for (size_t i = 0; i < chull.size()-1; ++i)
	{
		double dv = calc_triangle_area(chull[0], chull[i], chull[i + 1]);
		V += dv;
		res += dv*0.3333333333*Vector2D(chull[0].x + chull[i].x + chull[i + 1].x, chull[0].y + chull[i].y + chull[i + 1].y);
	}
	res *= 1.0/V;
	return res;
}

double SlabSymmetry::calcVolume(const vector<Vector2D>& chull) const
{
  double res = 0;
  const Vector2D anchor = chull[0];
  for(size_t i=1;i<chull.size()-1;++i)
    res += calc_triangle_area(anchor,chull[i],chull[i+1]);
  return res;
}

Axis::Axis(const Vector2D& origin_i,
	   const Vector2D& direction_i):
  origin(origin_i), direction(direction_i/abs(direction_i)) {}

namespace 
{
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
  return M_PI*(p1.y+p2.y)*sqrt((p2.y - p1.y )*(p2.y - p1.y )+ (p2.x - p1.x)*(p2.x - p1.x));
}

namespace {
  double calc_cone_segment_volume(const Vector2D& p1,
				  const Vector2D& p2,
				  const Axis& axis)
  {
    const Vector2D q1 = change_coordinate(p1, axis);
    const Vector2D q2 = change_coordinate(p2, axis);
    return (M_PI/3.)*(q2.x-q1.x)*(q1.y*q1.y+ q2.y *q2.y +q1.y*q2.y);
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

Vector2D CylindricalSymmetry::calcCentroid(const vector<Vector2D>& chull) const
{
  //	double V = 0;
	Vector2D res;
	Vector2D anchor = change_coordinate(chull[0], axis_);
	for (size_t i = 1; i < chull.size() - 1; ++i)
	{
		Vector2D q1 = change_coordinate(chull[i], axis_);
		Vector2D q2 = change_coordinate(chull[i+1], axis_);
		double dv = 0.5*std::abs(CrossProduct(q2 - anchor, q1 - anchor));
		//		V += dv;
		res += dv*Vector2D((2 * anchor.x*anchor.y + q1.x*anchor.y + q2.x*anchor.y + q1.y*anchor.x + q1.y*q1.x * 2 + q1.y*q2.x
			+ q2.y*anchor.x + q2.y*q1.x + 2 * q2.x*q2.y) / 12.0, (anchor.y*anchor.y + anchor.y*q1.y + anchor.y*q2.y
				+ q1.y*q2.y + q1.y*q1.y + q2.y*q2.y) / 6.0);
	}
	double realV = calcVolume(chull);
	res *= 2*M_PI / realV;
	
	Vector2D perp = zcross(axis_.direction);
	Vector2D res2(res.x*axis_.direction.x+res.y*perp.x, res.y*perp.y+res.x*axis_.direction.y);
	return res2;
}

double CylindricalSymmetry::calcVolume(const vector<Vector2D>& chull) const
{
  double res = 0;
  const Vector2D anchor = chull[0];
  for(size_t i=1;i<chull.size()-1;++i)
    res += calc_triangular_ring_volume(anchor,
				       chull[i],
				       chull[i+1],
				       axis_);
  return res;
}

const Axis& CylindricalSymmetry::getAxis(void) const
{
  return axis_;
}
