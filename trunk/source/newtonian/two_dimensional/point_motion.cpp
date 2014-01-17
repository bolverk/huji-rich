#include "point_motion.hpp"

vector<Vector2D> PointMotion::calcAllVelocities
(Tessellation const* tess,
 vector<Primitive> const& cells,
 double time)
{
  vector<Vector2D> res(cells.size());
  for(int i=0;i<(int)cells.size();++i)
    res[i] = CalcVelocity(i,tess,cells,time);
  return res;
}

PointMotion::~PointMotion(void) {}
