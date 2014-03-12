#include "edge_repellant.hpp"

EdgeRepellant::EdgeRepellant(PointMotion& naive,double inner_radius,
	double outer_radius,int total_specials,bool move_inner,int inner1,int inner2,
		int inner3,int n2start,int n3start,double inner_r1,double inner_r2,
	double inner_r3,double mass,double t0):
  naive_(naive),
  inner_radius_(inner_radius),
  outer_radius_(outer_radius),
  total_specials_(total_specials),
  move_inner_(move_inner),
  inner_n1_(inner1),
  inner_n2_(inner2),
  inner_n3_(inner3),
  n2start_(n2start),
  n3start_(n3start),
  inner_r1_(inner_r1),
  inner_r2_(inner_r2),
  inner_r3_(inner_r3),
  mass_(mass),t0_(t0),init_angles_(vector<double>()),first_time_(true)
  {}

Vector2D EdgeRepellant::CalcVelocity
(int index, 
 Tessellation const& tess,
 vector<Primitive> const& cells,
 double time)
{
  return naive_.CalcVelocity(index,tess,cells,time);
}

vector<Vector2D> EdgeRepellant::calcAllVelocities(Tessellation const& tess,
						  vector<Primitive> const& cells,
						  double time)
{
	if(move_inner_&&first_time_)
	  {
		  first_time_=false;
		  init_angles_.resize(total_specials_);
		  for(int i=0;i<total_specials_;++i)
		  {
			  const Vector2D point=tess.GetMeshPoint(i);
			  init_angles_[i]=atan2(point.y,point.x);
				  if(point.y<0)
					  init_angles_[i]+=2*M_PI;
		  }
	  }
  vector<Vector2D> result = naive_.calcAllVelocities
    (tess,cells,time);
  if(move_inner_)
  {
	  for(int i=0;i<total_specials_;++i)
	  {
		  const Vector2D point=tess.GetMeshPoint(i);
		  const double r=abs(point);
		  const Vector2D phi_hat(-point.y/r,point.x/r);
		  const Vector2D r_hat(point.x/r,point.y/r);
		  const double omega=sqrt(mass_/(inner_r2_*inner_r2_*inner_r2_));
		  const double vmag=omega*r;
		  double p_angle=atan2(point.y,point.x);
		  if(point.y<0)
			  p_angle+=2*M_PI;
		  double time_angle=(time-t0_)*omega+init_angles_[i];
		  time_angle-=2*M_PI*floor(time_angle/(2*M_PI));
		  double d_angle;
		  if(abs(time_angle-p_angle)<0.1)
				d_angle=time_angle-p_angle;
		  else
			  if(time_angle>p_angle)
				  d_angle=(time_angle-2*M_PI)-p_angle;
			  else
				  d_angle=time_angle-(p_angle-2*M_PI);
		  if(i<inner_n1_)
			result[i]=vmag*phi_hat+vmag*(1-r/inner_r1_)*r_hat+
			r*vmag*d_angle*phi_hat/tess.GetWidth(i);
		  else
			  if(i<n2start_)
				  result[i]=Vector2D(0,0);
			  else
				  if(i<(n2start_+inner_n2_))
					  result[i]=vmag*phi_hat+vmag*(1-r/inner_r2_)*r_hat+
					  r*vmag*d_angle*phi_hat/tess.GetWidth(i);
				  else
					  if(i<n3start_)
						  result[i]=Vector2D(0,0);
					  else
						  if(i<(n3start_+inner_n3_))
							  result[i]=vmag*phi_hat+vmag*(1-r/inner_r3_)*r_hat
							  +r*vmag*d_angle*phi_hat/tess.GetWidth(i);
						  else
							  result[i]=Vector2D(0,0);
	  }	  
  }
  else
  for(int i=0;i<total_specials_;++i)
    result[i] = Vector2D(0,0);
  for(int i=total_specials_;i<(int)result.size();++i)
  {
    const Vector2D mp = tess.GetMeshPoint(i);
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
