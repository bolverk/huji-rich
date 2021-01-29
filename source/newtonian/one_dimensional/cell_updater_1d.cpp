#include <cmath>
#include "cell_updater_1d.hpp"
#include "spdlog/spdlog.h"
#include "../../misc/serial_generate.hpp"

CellUpdater1D::~CellUpdater1D(void) {}

SimpleCellUpdater1D::SimpleCellUpdater1D(void) {}

namespace{

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
    res.tracers = serial_generate<double, double>
      (res.tracers,
       [&](double t){return t/extensive.mass;});
    return res;
  }

  void validate_input
    (const vector<Extensive>& extensives,
     const SimulationState1D& old)
  {
    assert(extensives.at(0).tracers.size()==
	   old.getCells().at(0).tracers.size());
  }

  template<class T> vector<T> diff(const vector<T> source){
    vector<T> res(source.size()-1);
    for(size_t i=0;i<res.size();++i)
      res.at(i) = source.at(i+1) - source.at(i);
    return res;
  }
}

vector<ComputationalCell> SimpleCellUpdater1D::operator()
  (const PhysicalGeometry1D& pg,
   const vector<Extensive>& extensives,
   const SimulationState1D& old,
   const EquationOfState& eos) const
{
  validate_input(extensives, old);
  const vector<double> volumes =
    diff(serial_generate<double,double>
	 (old.getVertices(),
	  [&](double r){return pg.calcVolume(r);}));
  vector<ComputationalCell> res(extensives.size());
  transform(extensives.begin(),
	    extensives.end(),
	    volumes.begin(),
	    res.begin(),
	    [&](const Extensive& e,
		const double& v){
	      return retrieve_single_cell(v,e,eos);});
  return res;
}

SimpleCellUpdater1D::~SimpleCellUpdater1D(void) {}
