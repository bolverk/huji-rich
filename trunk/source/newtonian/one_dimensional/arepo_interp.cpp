#include "arepo_interp.hpp"
#include "../common/hydrodynamics.hpp"
#include <cmath>

using namespace std;

ArepoInterp::ArepoInterp(IdealGas const& eos):
  eos_(eos) {}

namespace {
  Primitive Derivative
  (vector<double> const& grid,
   vector<Primitive> const& cells,
   size_t idx)
  {
    return (cells[size_t(idx)]-cells[size_t(idx-1)])/
      (0.5*(grid[size_t(idx+1)]-grid[size_t(idx-1)]));
  }

  Primitive SymmetricDerivative
  (vector<double> const& grid,
   vector<Primitive> const& cells,
   size_t idx)
  {
    const double xr = 0.5*(grid[size_t(idx+2)]+grid[size_t(idx+1)]);
    const double xl = 0.5*(grid[size_t(idx)]+grid[size_t(idx-1)]);
    return (cells[size_t(idx+1)]-cells[size_t(idx-1)])/(xr-xl);
  }

  bool effectively_zero(double x)
  {
    return abs(x)<1e-14;
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
				   size_t idx, int dir, double dt) const
{
  if(dir==0){
    if(idx==1)
      return hv[0];
    else if(idx==vp.size()-1)
      return hv[vp.size()-2];
    else if(idx==0)
      throw UniversalError("Error in PLM1D::InterpState\n Target of interpolation lies outside grid, on the left side");
    else{
      const Primitive slope_left = Derivative(vp,hv,idx-1);
      const Primitive slope_center = Derivative(vp,hv,idx);
      const Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx-1);
      const Primitive slope = minmod(slope_left,slope_center,2*slope_symmetric);
      const double dx = vp[size_t(idx)]-vp[size_t(idx-1)];
      const double edge_density = 
	hv[size_t(idx-1)].Density
	+slope.Density*dx/2
	-(dt/2)*slope.Density*(hv[size_t(idx-1)].Velocity.x-v_i)
	-(dt/2)*slope.Velocity.x*hv[size_t(idx-1)].Density;
      const double g = eos_.getAdiabaticIndex();
      const double edge_pressure = 
	hv[size_t(idx-1)].Pressure
	+slope.Pressure*dx/2
	-(dt/2)*slope.Velocity.x*g*hv[size_t(idx-1)].Pressure
	-(dt/2)*slope.Pressure*(hv[size_t(idx-1)].Velocity.x-v_i);
      const double edge_xvelocity = 
	hv[size_t(idx-1)].Velocity.x
	+(dx/2)*slope.Velocity.x
	-(dt/2)*slope.Velocity.x*(hv[size_t(idx-1)].Velocity.x-v_i)
	-(dt/2)*slope.Pressure/hv[size_t(idx-1)].Density;
      const double edge_yvelocity =
	hv[size_t(idx-1)].Velocity.y
	+(dx/2)*slope.Velocity.y
	-(dt/2)*slope.Velocity.y*(hv[size_t(idx-1)].Velocity.x-v_i);
      return CalcPrimitive
	(edge_density,
	 edge_pressure,
	 Vector2D(edge_xvelocity,
		  edge_yvelocity),
	 eos_);
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
      const Primitive slope_right = Derivative(vp,hv,idx+1);
      const Primitive slope_center = Derivative(vp,hv,idx);
      const Primitive slope_symmetric = SymmetricDerivative(vp,hv,idx);
      const Primitive slope = minmod(slope_center,slope_right,2*slope_symmetric);
      const double dx = vp[size_t(idx+1)]-vp[size_t(idx)];
      const double edge_density = 
	hv[size_t(idx)].Density
	-slope.Density*dx/2
	-(dt/2)*slope.Density*(hv[size_t(idx)].Velocity.x-v_i)
	-(dt/2)*slope.Velocity.x*hv[size_t(idx)].Density;
      const double g = eos_.getAdiabaticIndex();
      const double edge_pressure = 
	hv[size_t(idx)].Pressure
	-slope.Pressure*dx/2
	-(dt/2)*slope.Velocity.x*g*hv[size_t(idx)].Pressure
	-(dt/2)*slope.Pressure*(hv[size_t(idx)].Velocity.x-v_i);
      const double edge_xvelocity = 
	hv[size_t(idx)].Velocity.x
	-(dx/2)*slope.Velocity.x
	-(dt/2)*slope.Velocity.x*(hv[size_t(idx)].Velocity.x-v_i)
	-(dt/2)*slope.Pressure/hv[size_t(idx)].Density;
      const double edge_yvelocity =
	hv[size_t(idx)].Velocity.y
	-(dx/2)*slope.Velocity.y
	-(dt/2)*slope.Velocity.y*(hv[size_t(idx)].Velocity.x-v_i);
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
