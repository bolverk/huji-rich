#include "sb_edge_repellant.hpp"

SBEdgeRepellant::SBEdgeRepellant(PointMotion& naive,
				 double inner_radius1,
				 Vector2D const& center1,
				 double inner_radius2,
				 Vector2D const& center2,
				 double outer_radius,
				 int total_specials):
  naive_(naive),
  inner_radius1_(inner_radius1),
  center1_(center1),
  inner_radius2_(inner_radius2),
  center2_(center2),
  outer_radius_(outer_radius),
  total_specials_(total_specials) {}

Vector2D SBEdgeRepellant::CalcVelocity
(int index, 
 Tessellation const& tess,
 vector<Primitive> const& cells,
 double time)
{
  return naive_.CalcVelocity(index,tess,cells,time);
}

bool SBEdgeRepellant::should_be_decelerated(Vector2D const& mp,
					    Vector2D const& vm) const
{
  if(abs(mp)>outer_radius_)
    return ScalarProd(mp,vm)>0;
  else if(abs(mp-center1_)<inner_radius1_)
    return ScalarProd(mp-center1_,vm)<0;
  else if(abs(mp-center2_)<inner_radius2_)
    return ScalarProd(mp-center2_,vm)<0;
  else
    return false;
}

vector<Vector2D> SBEdgeRepellant::calcAllVelocities(Tessellation const& tess,
						  vector<Primitive> const& cells,
						  double time,vector<CustomEvolution*> &cevolve)
{
  vector<Vector2D> result = naive_.calcAllVelocities
    (tess,cells,time,cevolve);
  for(size_t i=0;i<static_cast<size_t>(total_specials_);++i)
    result[i] = Vector2D(0,0);
  for(size_t i=static_cast<size_t>(total_specials_);i<result.size();++i){
    const Vector2D mp = tess.GetMeshPoint(static_cast<int>(i));
    const Vector2D vm = cells[i].Velocity;
    if(should_be_decelerated(mp,vm)){
      const Vector2D rhat = mp/abs(mp);
      const double vr = ScalarProd(rhat,vm);
      result[i] = 0.1*vr*rhat+(vm-vr*rhat);
    }
  }
  return result;
}
