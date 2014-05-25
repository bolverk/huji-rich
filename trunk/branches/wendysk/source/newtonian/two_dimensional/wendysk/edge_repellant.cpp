#include "edge_repellant.hpp"

EdgeRepellant::EdgeRepellant(PointMotion& naive,
			     double inner_radius,
			     double outer_radius,
			     int total_specials):
  naive_(naive),
  inner_radius_(inner_radius),
  outer_radius_(outer_radius),
  total_specials_(total_specials){}

Vector2D EdgeRepellant::CalcVelocity
(int index, 
 Tessellation const* tess,
 vector<Primitive> const& cells,
 double time)
{
  return naive_.CalcVelocity(index,tess,cells,time);
}

vector<Vector2D> EdgeRepellant::calcAllVelocities(Tessellation const* tess,
						  vector<Primitive> const& cells,
						  double time)
{
  vector<Vector2D> result = naive_.calcAllVelocities
    (tess,cells,time);
  for(int i=0;i<total_specials_;++i)
    result[i] = Vector2D(0,0);
  for(int i=total_specials_;i<(int)result.size();++i){
    const Vector2D mp = tess->GetMeshPoint(i);
    const Vector2D vm = cells[i].Velocity;
    if((ScalarProd(mp,vm)>0&&abs(mp)>outer_radius_)||
       (ScalarProd(mp,vm)<0&&abs(mp)<inner_radius_)){
      const Vector2D rhat = mp/abs(mp);
      const double vr = ScalarProd(rhat,vm);
      result[i] = 0.1*vr*rhat+(vm-vr*rhat);
    }
  }
  return result;
}
