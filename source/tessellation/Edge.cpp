#include "Edge.hpp"
#include <cmath>
#include "../misc/universal_error.hpp"

Edge::Edge(void):
vertices(Vector2D(), Vector2D()),
	neighbors(0,0) {}

Edge& Edge::operator=(const Edge& other)
{
  vertices = other.vertices;
  neighbors = other.neighbors;
  return *this;
}

Edge::~Edge(void){}

Edge::Edge(Edge const& other):
vertices(other.vertices),
	neighbors(other.neighbors) {}

Edge::Edge(Vector2D const& p1, Vector2D const& p2,
	int neighbor1, int neighbor2):
vertices(p1,p2), neighbors(neighbor1, neighbor2) {}

double Edge::GetLength(void) const
{
	return abs(vertices.second-vertices.first);
}

double DistanceToEdge(Vector2D const& point,Edge const& edge)
{
	Vector2D v=edge.vertices.second-edge.vertices.first;
	Vector2D w=point-edge.vertices.first;
	double c1,c2;
	c1=ScalarProd(v,w);
	if(c1<=0)
		return point.distance(edge.vertices.first);
	c2=ScalarProd(v,v);
	if(c2<=c1)
		return point.distance(edge.vertices.second);
	return point.distance(edge.vertices.first+(c1/c2)*v);
}

bool SegmentIntersection(Edge const&edge1,Edge const&edge2,
	Vector2D &Intersection,double eps)
{
	bool res=true;
	const double areascale=std::min((edge1.vertices.first.x-edge1.vertices.second.x)*
		(edge1.vertices.first.x-edge1.vertices.second.x)+(edge1.vertices.first.y-
		edge1.vertices.second.y)*(edge1.vertices.first.y-edge1.vertices.second.y),
		(edge2.vertices.first.x-edge2.vertices.second.x)*
		(edge2.vertices.first.x-edge2.vertices.second.x)+(edge2.vertices.first.y-
		edge2.vertices.second.y)*(edge2.vertices.first.y-edge2.vertices.second.y));
	if(std::min(edge1.vertices.second.x,edge1.vertices.first.x)>std::max(edge2.vertices.second.x,edge2.vertices.first.x)||
	   std::min(edge2.vertices.second.x,edge2.vertices.first.x)>std::max(edge1.vertices.second.x,edge1.vertices.first.x)||
	   std::min(edge1.vertices.second.y,edge1.vertices.first.y)>std::max(edge2.vertices.second.y,edge2.vertices.first.y)||
	   std::min(edge2.vertices.second.y,edge2.vertices.first.y)>std::max(edge1.vertices.second.y,edge1.vertices.first.y))
		res=false;
	double d=(edge1.vertices.first.x-edge1.vertices.second.x)*(edge2.vertices.first.y-edge2.vertices.second.y)
		-(edge2.vertices.first.x-edge2.vertices.second.x)*(edge1.vertices.first.y-edge1.vertices.second.y);
	if(fabs(d)<1e-8*areascale)
		return false;
	double xi=((edge2.vertices.first.x-edge2.vertices.second.x)*(edge1.vertices.first.x*edge1.vertices.second.y-
		edge1.vertices.second.x*edge1.vertices.first.y)-(edge1.vertices.first.x-edge1.vertices.second.x)*
		(edge2.vertices.first.x*edge2.vertices.second.y-edge2.vertices.second.x*edge2.vertices.first.y))/d;
	double yi=((edge2.vertices.first.y-edge2.vertices.second.y)*(edge1.vertices.first.x*edge1.vertices.second.y-
		edge1.vertices.second.x*edge1.vertices.first.y)-(edge1.vertices.first.y-edge1.vertices.second.y)*
		(edge2.vertices.first.x*edge2.vertices.second.y-edge2.vertices.second.x*edge2.vertices.first.y))/d;
	Intersection.Set(xi,yi);
	eps=eps*sqrt(areascale);
	if((xi+eps)<std::min(edge1.vertices.first.x,edge1.vertices.second.x)||(xi-eps)>std::max(edge1.vertices.first.x,edge1.vertices.second.x))
		return false;
	if((xi+eps)<std::min(edge2.vertices.first.x,edge2.vertices.second.x)||(xi-eps)>std::max(edge2.vertices.first.x,edge2.vertices.second.x))
		return false;
	if((yi+eps)<std::min(edge1.vertices.first.y,edge1.vertices.second.y)||(yi-eps)>std::max(edge1.vertices.first.y,edge1.vertices.second.y))
		return false;
	if((yi+eps)<std::min(edge2.vertices.first.y,edge2.vertices.second.y)||(yi-eps)>std::max(edge2.vertices.first.y,edge2.vertices.second.y))
		return false;
	return res;
}

Vector2D Parallel(Edge const& edge)
{
	return (edge.vertices.second - edge.vertices.first);
}

Vector2D calc_centroid(const Edge& edge)
{
  return 0.5*(edge.vertices.first+edge.vertices.second);
}
