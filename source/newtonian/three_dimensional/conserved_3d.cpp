#include "conserved_3d.hpp"

Conserved3D::Conserved3D(void):
  mass(0), momentum(), energy(0), tracers() {}

Conserved3D::Conserved3D(double mass_i,
			 const Vector3D& momentum_i,
			 double energy_i):
  mass(mass_i), momentum(momentum_i), energy(energy_i), tracers() {}

Conserved3D::Conserved3D(double mass_i,
			 const Vector3D& momentum_i,
			 double energy_i,
			 const vector<double>& tracers_i):
  mass(mass_i), momentum(momentum_i), 
  energy(energy_i), tracers(tracers_i) {}

namespace {
  vector<double> operator*(double s, const vector<double>& v)
  {
    vector<double> res(v.size());
    for(size_t i=0;i<v.size();++i)
      res[i] = s*v[i];
    return res;
  }
}

Conserved3D operator*(double s, const Conserved3D& c)
{
  return Conserved3D(s*c.mass,
		     s*c.momentum,
		     s*c.energy,
		     s*c.tracers);
}
