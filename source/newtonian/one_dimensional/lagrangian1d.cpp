#include "lagrangian1d.hpp"

Lagrangian1D::Lagrangian1D(bool rigid_walls):
  rigid_walls_(rigid_walls) {}

double Lagrangian1D::operator()
(size_t i, vector<double> const& /*vp*/,
 vector<ComputationalCell> const& hv) const
{
  if(i==0){
    if(rigid_walls_)
      return 0;
    else
      return hv[0].velocity.x;
  }
  else if(i==hv.size()){
    if(rigid_walls_)
      return 0;
    else
      return hv[hv.size()-1].velocity.x;
  }
  else
    return 0.5*(hv[size_t(i)-1].velocity.x+
		hv[size_t(i)].velocity.x);
}
