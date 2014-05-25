#include "eos_consistent1d.hpp"

using namespace std;
using namespace interpolations1d;

EOSConsistent::EOSConsistent
(SpatialReconstruction1D const& naive,
 EquationOfState const& eos):
  naive_(naive), eos_(eos) {}

Primitive EOSConsistent::InterpState
(vector<double> const& vp,
 vector<Primitive> const& hv,
 double interface_speed,
 size_t i, int dir, double dt) const
{
  Primitive res = naive_.InterpState
    (vp,hv,interface_speed,i,dir,dt);
  res.Energy = eos_.dp2e(res.Density,res.Pressure);
  res.SoundSpeed = eos_.dp2c(res.Density,res.Pressure);
  return res;
}
