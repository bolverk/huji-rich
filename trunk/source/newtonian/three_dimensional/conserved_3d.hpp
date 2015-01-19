#ifndef CONSERVED_3D_HPP
#define CONSERVED_3D_HPP 1

#include "../../3D/GeometryCommon/Vector3D.hpp"

//! \brief Conserved variables for a 3D computational cell
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

  Conserved3D& operator-=(const Conserved3D& diff);

  Conserved3D& operator+=(const Conserved3D& diff);
};

Conserved3D operator*(double s, const Conserved3D& c);

Conserved3D operator/(const Conserved3D& c, double s);

#endif // CONSERVED_3D_HPP
