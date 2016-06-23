#include "computational_cell.hpp"

ComputationalCell3D::ComputationalCell3D(void):
  density(0), pressure(0), velocity(), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i):
  density(density_i), pressure(pressure_i), 
  velocity(velocity_i), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i,
				     const vector<double>& tracers_i,
					 const vector<bool>& stickers_i):
  density(density_i), pressure(pressure_i),
  velocity(velocity_i), tracers(tracers_i),stickers(stickers_i) {}

Slope3D::Slope3D(void) : xderivative(ComputationalCell3D()), yderivative(ComputationalCell3D()) {}

Slope3D::Slope3D(ComputationalCell3D const & x, ComputationalCell3D const & y) : xderivative(x), yderivative(y)
{}
