#include "point_motion.hpp"

vector<Vector2D> PointMotion::calcAllVelocities
(Tessellation const& tess,
 vector<Primitive> const& cells,
 double time,vector<CustomEvolution*> &cevolve)
{
  vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
  for(size_t i=0;i<static_cast<size_t>(tess.GetPointNo());++i)
	  if(cevolve[i]==0)
	    res[i] = CalcVelocity(static_cast<int>(i),tess,cells,time);
	  else
	    res[i]=cevolve[i]->CalcVelocity(static_cast<int>(i),tess,cells,time);
  return res;
}

PointMotion::~PointMotion(void) {}
