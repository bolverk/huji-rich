#include "Edge.hpp"

#include <cmath>

Edge::Edge(void): 
  _p1(Vector2D()),
  _p2(Vector2D()),
  _neighbor1(0), 
  _neighbor2(0) {}

Edge::~Edge(void){}

Edge::Edge(Edge const& other):
  _p1(other._p1), _p2(other._p2),
  _neighbor1(other._neighbor1),
  _neighbor2(other._neighbor2){}

Edge::Edge(Vector2D const& p1, Vector2D const& p2,
	   int neighbor1, int neighbor2):
  _p1(p1), _p2(p2), _neighbor1(neighbor1), _neighbor2(neighbor2) {}

int Edge::GetNeighbor(int index) const
{
  if(index==0)
    return _neighbor1;
  else if(index==1)
    return _neighbor2;
  else{
    InvalidSideIndexException eo(index);
    throw eo;
  }
}

double Edge::get_x(int index) const
{
  if(index==0) 
    return _p1.get_x();
  else 
    return _p2.get_x();
}

double Edge::get_y(int index) const
{
  if(index==0) 
    return _p1.get_y();
  else 
    return _p2.get_y();
}

/*int Edge::get_other_adjacent_cell(int cell) const
  {
  if(cell == _neighbor1)
  return _neighbor2;
  else if(cell == _neighbor2)
  return _neighbor1;
  else
  {
  InvalidCellEdgeException eo(cell);
  throw eo;
  }
  }*/


Vector2D Edge::GetVertex(int index) const
{
  if(index==0)
    return _p1;
  else if(index==1)
    return _p2;
  else{
    InvalidSideIndexException eo(index);
    throw eo;
  }
}

double Edge::GetLength(void) const
{
  return abs(_p1-_p2);
}

void Edge::set_friend(int dim,int data)
{
  if(dim==0)
    _neighbor1=data;
  else
    _neighbor2=data;
}

void Edge::set_x(int p,double data)
{
  if(p==0)
    _p1.x=data;
  else
    _p2.x=data;
}

void Edge::set_y(int p,double data)
{
  if(p==0)
    _p1.y=data;
  else
    _p2.y=data;
}

InvalidSideIndexException::InvalidSideIndexException(int index)
  :	_index(index)
{}


/*int InvalidSideIndexException::GetIndex(void) const
  {
  return _index;
  }*/


InvalidCellEdgeException::InvalidCellEdgeException(int index)
  :	_index(index)
{}

/*int InvalidCellEdgeException::get_index(void) const
  {
  return _index;
  }*/

double DistanceToEdge(Vector2D const& point,Edge const& edge)
{
  Vector2D v=edge.GetVertex(1)-edge.GetVertex(0);
  Vector2D w=point-edge.GetVertex(0);
  double c1,c2;
  c1=ScalarProd(v,w);
  if(c1<=0)
    return point.distance(edge.GetVertex(0));
  c2=ScalarProd(v,v);
  if(c2<=c1)
    return point.distance(edge.GetVertex(1));
  return point.distance(edge.GetVertex(0)+(c1/c2)*v);
}

bool SegmentIntersection(Edge const& edge1,Edge const& edge2,
			 Vector2D &Intersection)
{
  double d=(edge1.get_x(0)-edge1.get_x(1))*(edge2.get_y(0)-edge2.get_y(1))
    -(edge2.get_x(0)-edge2.get_x(1))*(edge1.get_y(0)-edge1.get_y(1));
  if(d==0)
    return false;
  double xi=((edge2.get_x(0)-edge2.get_x(1))*(edge1.get_x(0)*edge1.get_y(1)-
					      edge1.get_x(1)*edge1.get_y(0))-(edge1.get_x(0)-edge1.get_x(1))*
	     (edge2.get_x(0)*edge2.get_y(1)-edge2.get_x(1)*edge2.get_y(0)))/d;
  double yi=((edge2.get_y(0)-edge2.get_y(1))*(edge1.get_x(0)*edge1.get_y(1)-
					      edge1.get_x(1)*edge1.get_y(0))-(edge1.get_y(0)-edge1.get_y(1))*
	     (edge2.get_x(0)*edge2.get_y(1)-edge2.get_x(1)*edge2.get_y(0)))/d;
  Intersection.Set(xi,yi);
  double eps=1e-8*min(edge1.GetLength(),edge2.GetLength());
  if((xi+eps)<min(edge1.get_x(0),edge1.get_x(1))||(xi-eps)>max(edge1.get_x(0),edge1.get_x(1)))
    return false;
  if((xi+eps)<min(edge2.get_x(0),edge2.get_x(1))||(xi-eps)>max(edge2.get_x(0),edge2.get_x(1)))
    return false;
  return true;
}

Vector2D Parallel(Edge const& edge)
{
  return (edge.GetVertex(1) - edge.GetVertex(0));
}

