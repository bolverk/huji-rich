#include "triangle.hpp"
// Implementation based on 
// http://www.blackpawn.com/texts/pointinpoly/default.html

namespace {
bool SameSide(Vector2D p1, Vector2D p2, Vector2D a, Vector2D b)
{
  double cp1 = CrossProduct(b-a,p1-a);
  double cp2 = CrossProduct(b-a,p2-a);
  return cp1*cp2>=0;
}
}

Triangle::Triangle(vector<Vector2D> vv, double vi, double vo):
  vv_(vv), vi_(vi), vo_(vo) {}

double Triangle::EvalAt(Vector2D const& r) const
{
  if (SameSide(r,vv_[0],vv_[1],vv_[2]) &&
    SameSide(r,vv_[1],vv_[0],vv_[2]) &&
      SameSide(r,vv_[2],vv_[0],vv_[1]))
    return vi_;
  else
    return vo_;
}
