#ifndef CONSERVED_3D_HPP
#define CONSERVED_3D_HPP 1

#include "../../3D/Tessellation/Vector3D.hpp"

class Conserved3D
{
public:

  double mass;
  Vector3D momentum;
  double energy;
  vector<double> tracers;

  Conserved3D(void);

  Conserved3D(double mass_i,
	      const Vector3D& momentum_i,
	      double energy_i);

  Conserved3D(double mass_i,
	      const Vector3D& momentum_i,
	      double energy_i,
	      const vector<double>& tracers_i);
};

#endif // CONSERVED_3D_HPP
