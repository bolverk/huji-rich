#include "pcm1d.hpp"
#include "../../misc/universal_error.hpp"

Primitive PCM1D::InterpState(vector<double> const& /*vp*/,
			     vector<Primitive> const& hv,
			     double /*interface_speed*/,
			     size_t i, int dir, double /*dt*/) const
{
  if(dir==0)//left boundary 
    return hv[i-1];
  else if(dir==1) // right boundary
    return hv[i];
  else
    throw UniversalError("Invalid value for argument dir in PCM1D::InterpState");
}
