#include "point_motion.hpp"

vector<Vector2D> PointMotion::calcAllVelocities
(Tessellation const& tess,
 vector<Primitive> const& cells,
 double time,vector<CustomEvolution*> &cevolve)
{
  vector<Vector2D> res((size_t)tess.GetPointNo());
  for(int i=0;i<tess.GetPointNo();++i)
	  if(cevolve[i]==0)
		res[i] = CalcVelocity(i,tess,cells,time);
	  else
		  res[i]=cevolve[i]->CalcVelocity(i,tess,cells,time);
  return res;
}

PointMotion::~PointMotion(void) {}
