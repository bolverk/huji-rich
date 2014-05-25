#include <cmath>
#include "plm1d.hpp"
#include "../../misc/universal_error.hpp"

using namespace std;

namespace {
  Primitive Derivative(vector<double> const& grid,
		       vector<Primitive> const& cells,
		       size_t idx)
  {
    return (cells[idx]-cells[idx-1])/
      (0.5*(grid[idx+1]-grid[idx-1]));

  }

  Primitive SymmetricDerivative(vector<double> const& grid,
				vector<Primitive> const& cells,
				size_t idx)
  {
    const double xr = 0.5*(grid[idx+2]+grid[idx+1]);
    const double xl = 0.5*(grid[idx]+grid[idx-1]);
    return (cells[idx+1]-cells[idx-1])/(xr-xl);
  }

  bool effectively_zero(double x)
  {
    return fabs(x)<1e-14;
  }

  double sgn(double x)
  {
    if(x>0)
      return 1;
    else if(x<0)
      return -1;
    else if(effectively_zero(x))
      return 0;
    else
      throw UniversalError("NaN detected in sgn");
  }
}

namespace {
  double my_min(double x, double y, double z)
  {
    return min(min(x,y),z);
  }

  double minmod(double x, double y, double z)
  {
    return 0.125*(sgn(x)+sgn(y))*
      (sgn(y)+sgn(z))*
      (sgn(z)+sgn(x))*
      my_min(abs(x),abs(y),abs(z));
  }

  Primitive minmod(Primitive const& p1,
		   Primitive const& p2,
		   Primitive const& p3)
  {
    Primitive res;
    for(int i=0;i<res.GetVarNo();++i)
      res[i] = minmod(p1[i],p2[i],p3[i]);
    return res;
  }
}

PLM1D::PLM1D(bool second_order_time,
	     bool slope_limiter_flag):
  second_order_time_(second_order_time),
  slope_limiter_flag_(slope_limiter_flag) {}

Primitive PLM1D::InterpState(vector<double> const& vp,
			     vector<Primitive> const& hv,
			     double /*interface_speed*/,
			     size_t idx, 
			     int dir, double dt) const
{
  if(dir==0){
    if(idx==1)
      return hv[0];
    else if(idx==vp.size()-1)
      return hv[vp.size()-2];
    else if(idx==0)
      throw UniversalError("Error in PLM1D::InterpState\n Target of interpolation lies outside grid, on the left side");
    else{
      Primitive slope_left = Derivative(vp,hv,idx-1);
      Primitive slope_center = Derivative(vp,hv,idx);
      Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx-1);
      Primitive slope = minmod(slope_left,slope_center,2*slope_symmetric);
      double dx = vp[idx]-vp[idx-1];
      if(!slope_limiter_flag_)
	return hv[idx-1]+slope_center*dx/2;
      if(second_order_time_)
	return hv[idx-1]+slope*(dx/2-hv[idx-1].SoundSpeed*dt/2);
      else
	return hv[idx-1]+slope*dx/2;
    }
  }
  else if(dir==1){
    if(idx==vp.size()-2)
      return hv[vp.size()-2];
    else if(idx==0)
      return hv[0];
    else if(idx==vp.size()-1)
      throw UniversalError("Error in PLM1D::InterpState\n Target of interpolation lies outside grid, on the right side");
    else{
      Primitive slope_right = Derivative(vp,hv,idx+1);
      Primitive slope_center = Derivative(vp,hv,idx);
      Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx);
      Primitive slope = minmod(slope_center,slope_right,2*slope_symmetric);
      double dx = vp[idx+1]-vp[idx];
      if(!slope_limiter_flag_)
	return hv[idx]-slope_center*dx/2;
      if(second_order_time_)
	return hv[idx]-slope*(dx/2-hv[idx].SoundSpeed*dt/2);
      else
	return hv[idx]-slope*dx/2;
    }
  }
  else
    throw UniversalError("Unknown direction in plm1d");
}
