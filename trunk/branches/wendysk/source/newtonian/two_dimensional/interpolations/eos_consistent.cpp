#include "eos_consistent.hpp"

EOSConsistent::EOSConsistent(SpatialReconstruction& sr,
			     EquationOfState const& eos):
  sr_(sr), eos_(eos) {}

void EOSConsistent::Prepare(Tessellation const* tess,vector<Primitive> 
	const& cells,double dt,vector<bool> const& mask,double time)
{
  sr_.Prepare(tess,cells,dt,mask,time);
}

Primitive EOSConsistent::Interpolate(Tessellation const* tess,
	vector<Primitive> const& cells,double dt,Edge const& edge,int side,
	InterpolationType interptype) const
{
  Primitive res = sr_.Interpolate(tess,cells,dt,edge,side,interptype);
  res.Energy = eos_.dp2e(res.Density,res.Pressure);
  res.SoundSpeed = eos_.dp2c(res.Density,res.Pressure);
  return res;
}

bool EOSConsistent::WasSlopeLimited(int index)const
{
	return sr_.WasSlopeLimited(index);
}
