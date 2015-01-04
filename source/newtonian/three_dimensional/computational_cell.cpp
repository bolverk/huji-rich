#include "computational_cell.hpp"

ComputationalCell::ComputationalCell(void):
  density(0), pressure(0), velocity(), tracers() {}

ComputationalCell::ComputationalCell(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i):
  density(density_i), pressure(pressure_i), 
  velocity(velocity_i), tracers() {}

ComputationalCell::ComputationalCell(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i,
				     const vector<double>& tracers_i):
  density(density_i), pressure(pressure_i),
  velocity(velocity_i), tracers(tracers_i) {}
