#include "lagrangian1d.hpp"

Lagrangian1D::Lagrangian1D(bool rigid_walls):
  rigid_walls_(rigid_walls) {}

double Lagrangian1D::CalcVelocity
(int i, vector<double> const& /*vp*/,
 vector<Primitive> const& hv) const
{
  if(i==0){
    if(rigid_walls_)
      return 0;
    else
      return hv[0].Velocity.x;
  }
  else if(i==static_cast<int>(hv.size())){
    if(rigid_walls_)
      return 0;
    else
      return hv[hv.size()-1].Velocity.x;
  }
  else
    return 0.5*(hv[size_t(i)-1].Velocity.x+
		hv[size_t(i)].Velocity.x);
}
