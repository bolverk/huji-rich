#include "arepo_interp.hpp"
#include "../common/hydrodynamics.hpp"

ArepoInterp::ArepoInterp(IdealGas const& eos):
  eos_(eos) {}

namespace {
  Primitive Derivative
  (vector<double> const& grid,
   vector<Primitive> const& cells,
   int idx)
  {
    return (cells[idx]-cells[idx-1])/
      (0.5*(grid[idx+1]-grid[idx-1]));
  }

  Primitive SymmetricDerivative
  (vector<double> const& grid,
   vector<Primitive> const& cells,
   int idx)
  {
    const double xr = 0.5*(grid[idx+2]+grid[idx+1]);
    const double xl = 0.5*(grid[idx]+grid[idx-1]);
    return (cells[idx+1]-cells[idx-1])/(xr-xl);
  }

  double sgn(double x)
  {
    if(x>0)
      return 1;
    else if(x<0)
      return -1;
    else if(x==0)
      return 0;
    else
      throw UniversalError("NaN detected in sgn");
  }

  double minmod(double x, double y)
  {
    return 0.5*(sgn(x)+sgn(y))*
      min(abs(x),abs(y));
  }

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

Primitive ArepoInterp::InterpState(vector<double> const& vp,
				   vector<Primitive> const& hv,
				   double v_i,
				   int idx, int dir, double dt) const
{
  if(dir==0){
    if(idx==1)
      return hv[0];
    else if(idx==(int)vp.size()-1)
      return hv[vp.size()-2];
    else if(idx==0)
      throw UniversalError("Error in PLM1D::InterpState\n Target of interpolation lies outside grid, on the left side");
    else{
      const Primitive slope_left = Derivative(vp,hv,idx-1);
      const Primitive slope_center = Derivative(vp,hv,idx);
      const Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx-1);
      const Primitive slope = minmod(slope_left,slope_center,2*slope_symmetric);
      const double dx = vp[idx]-vp[idx-1];
      const double edge_density = 
	hv[idx-1].Density
	+slope.Density*dx/2
	-(dt/2)*slope.Density*(hv[idx-1].Velocity.x-v_i)
	-(dt/2)*slope.Velocity.x*hv[idx-1].Density;
      const double g = eos_.getAdiabaticIndex();
      const double edge_pressure = 
	hv[idx-1].Pressure
	+slope.Pressure*dx/2
	-(dt/2)*slope.Velocity.x*g*hv[idx-1].Pressure
	-(dt/2)*slope.Pressure*(hv[idx-1].Velocity.x-v_i);
      const double edge_xvelocity = 
	hv[idx-1].Velocity.x
	+(dx/2)*slope.Velocity.x
	-(dt/2)*slope.Velocity.x*(hv[idx-1].Velocity.x-v_i)
	-(dt/2)*slope.Pressure/hv[idx-1].Density;
      const double edge_yvelocity =
	hv[idx-1].Velocity.y
	+(dx/2)*slope.Velocity.y
	-(dt/2)*slope.Velocity.y*(hv[idx-1].Velocity.x-v_i);
      return CalcPrimitive
	(edge_density,
	 edge_pressure,
	 Vector2D(edge_xvelocity,
		  edge_yvelocity),
	 eos_);
    }
  }
  else if(dir==1){
    if(idx==(int)vp.size()-2)
      return hv[vp.size()-2];
    else if(idx==0)
      return hv[0];
    else if(idx==(int)vp.size()-1)
      throw UniversalError("Error in PLM1D::InterpState\n Target of interpolation lies outside grid, on the right side");
    else{
      const Primitive slope_right = Derivative(vp,hv,idx+1);
      const Primitive slope_center = Derivative(vp,hv,idx);
      const Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx);
      const Primitive slope = minmod(slope_center,slope_right,2*slope_symmetric);
      const double dx = vp[idx+1]-vp[idx];
      const double edge_density = 
	hv[idx].Density
	-slope.Density*dx/2
	-(dt/2)*slope.Density*(hv[idx].Velocity.x-v_i)
	-(dt/2)*slope.Velocity.x*hv[idx].Density;
      const double g = eos_.getAdiabaticIndex();
      const double edge_pressure = 
	hv[idx].Pressure
	-slope.Pressure*dx/2
	-(dt/2)*slope.Velocity.x*g*hv[idx].Pressure
	-(dt/2)*slope.Pressure*(hv[idx].Velocity.x-v_i);
      const double edge_xvelocity = 
	hv[idx].Velocity.x
	-(dx/2)*slope.Velocity.x
	-(dt/2)*slope.Velocity.x*(hv[idx].Velocity.x-v_i)
	-(dt/2)*slope.Pressure/hv[idx].Density;
      const double edge_yvelocity =
	hv[idx].Velocity.y
	-(dx/2)*slope.Velocity.y
	-(dt/2)*slope.Velocity.y*(hv[idx].Velocity.x-v_i);
      return CalcPrimitive
	(edge_density,
	 edge_pressure,
	 Vector2D(edge_xvelocity,
		  edge_yvelocity),
	 eos_);
    }
  }
  else
    throw UniversalError("Unknown direction in plm1d");
}
