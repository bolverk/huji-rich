#include "default_cell_updater.hpp"

DefaultCellUpdater::DefaultCellUpdater(void) {}

ComputationalCell DefaultCellUpdater::operator()
(const Conserved3D& intensive,
 const EquationOfState& eos) const
{
  const double density = intensive.mass;
  const Vector3D velocity = intensive.momentum/density;
  const double thermal_energy = 
    intensive.energy - 0.5*density*ScalarProd(velocity,velocity);
  const double pressure = eos.de2p(density, thermal_energy/density);
  return ComputationalCell(density,pressure,velocity,intensive.tracers);
}
