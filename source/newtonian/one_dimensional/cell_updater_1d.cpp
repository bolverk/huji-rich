#include <cmath>
#include "cell_updater_1d.hpp"

CellUpdater1D::~CellUpdater1D(void) {}

SimpleCellUpdater1D::SimpleCellUpdater1D(void) {}

namespace{

  double calc_cell_volume
    (const PhysicalGeometry1D& pg,
     const vector<double>& vertices,
     const size_t i)
  {
    return pg.calcVolume(vertices.at(i+1)) - pg.calcVolume(vertices.at(i));
  }

  ComputationalCell retrieve_single_cell
    (const double volume,
     const Extensive& extensive,
     const EquationOfState& eos)
  {
    const double density = extensive.mass/volume;
    const Vector2D velocity = extensive.momentum/extensive.mass;
    const double kinetic_specific_energy = 0.5*pow(abs(velocity),2);
    const double total_specific_energy = extensive.energy/extensive.mass;
    const double thermal_specific_energy =
      total_specific_energy - kinetic_specific_energy;
    const double pressure = eos.de2p
      (density, thermal_specific_energy);

    ComputationalCell res;
    res.density = density;
    res.pressure = pressure;
    res.velocity = velocity;
    res.tracers = extensive.tracers;
    for(size_t i=0;i<res.tracers.size();++i)
      res.tracers.at(i) /= extensive.mass;
    return res;
  }

  /*
  ComputationalCell retrieve_single_cell
    (const Conserved& intensive,
     const EquationOfState& eos)
     {
         const double density = intensive.Mass;
         const Vector2D velocity = intensive.Momentum/intensive.Mass;
         const double kinetic_energy = 0.5*pow(abs(velocity),2);
         const double total_energy = intensive.Energy/intensive.Mass;
         const double thermal_energy = total_energy - kinetic_energy;
         const double pressure = eos.de2p(density, thermal_energy);
	 //         const double sound_speed = eos.dp2c(density, pressure);
	 ComputationalCell res;
	 res.density = density;
	 res.pressure = pressure;
	 res.velocity = velocity;
	 return res;
     }
  */
}

vector<ComputationalCell> SimpleCellUpdater1D::operator()
  (const PhysicalGeometry1D& pg,
   const vector<Extensive>& extensives,
   const SimulationState1D& old,
   const EquationOfState& eos) const
{
    vector<ComputationalCell> res(extensives.size());
    for(size_t i=0;i<res.size();++i)
        res.at(i) = retrieve_single_cell
	  (calc_cell_volume(pg,old.getVertices(),i),
	   extensives.at(i),
	   eos);
    return res;
}

SimpleCellUpdater1D::~SimpleCellUpdater1D(void) {}
