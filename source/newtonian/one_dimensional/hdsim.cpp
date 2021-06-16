#include <cmath>
#include <algorithm>
#include "hdsim.hpp"
#include "../../misc/universal_error.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#include "spdlog/spdlog.h"
#include "../../misc/serial_generate.hpp"

using namespace std;

// Diagnostics

double hdsim1D::GetTime(void) const
{
  return time_;
}

int hdsim1D::GetCycle(void) const
{
  return cycle_;
}

const SimulationState1D& hdsim1D::getState(void) const
{
  return ss_;
}

SimulationState1D& hdsim1D::getState(void)
{
  return ss_;
}

const vector<Extensive>& hdsim1D::getExtensives(void) const
{
  return extensives_;
}

void hdsim1D::setExtensives(const vector<Extensive>& new_ext)
{
  extensives_ = new_ext;
}

// External functions

namespace {

  Extensive calc_single_extensive
  (const ComputationalCell& cell,
   const double& volume,
   const EquationOfState& eos)
  {
    Extensive res;
    res.mass = cell.density*volume;
    res.momentum = res.mass*cell.velocity;
    const double kinetic_specific_energy =
      0.5*pow(abs(cell.velocity),2);
    const double thermal_specific_energy =
      eos.dp2e(cell.density, cell.pressure);
    res.energy = res.mass*
      (kinetic_specific_energy+thermal_specific_energy);
    res.tracers = serial_generate<double, double>
      (cell.tracers,
       [&](const double& t){return t*res.mass;});
    return res;
  }

  vector<Extensive> calc_extensives
  (const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   const EquationOfState& eos)
  {
    const vector<double>& vertices = ss.getVertices();
    const vector<ComputationalCell>& cells = ss.getCells();
    const vector<double> volumes = diff
      (serial_generate<double, double>
       (vertices, [&](const double& r){return pg.calcVolume(r);}));
    return serial_generate<ComputationalCell, double, Extensive>
      (cells,
       volumes,
       [&](const ComputationalCell& c,
	   const double& v){
	 return calc_single_extensive(c,v,eos);});
  }

#if 0
// ***Start*** Added by Emma	
ComputationalCell calc_single_simulationstate
  (const Extensive& ext,
   const double& volume,
   const EquationOfState& eos)
  {
    SimulationState1D cell;
    cell.density = ext.mass/volume;
    cell.velocity = ext.momentum/ext.mass;
    const double kinetic_specific_energy = 0.5*pow(abs(cell.velocity), 2);
    const double thermal_specific_energy = ext.energy/ext.mass - kinetic_specific_energy;
    cell.pressure = eos.de2p(cell.density, thermal_specific_energy);
    return cell;

  }

  SimulationState1D calc_simulationstates
  (const PhysicalGeometry1D& pg,
   const vector<Extensives>& ext,
   const EquationOfState& eos)
  {
    const vector<double>& vertices = ss.getVertices();
    const vector<ComputationalCell>& cells = ss.getCells();
    const vector<double> volumes = diff
      (serial_generate<double, double>
       (vertices, [&](const double& r){return pg.calcVolume(r);}));
    return serial_generate<ComputationalCell, double, Extensive>
      (cells,
       volumes,
       [&](const ComputationalCell& c,
	   const double& v){
	 return calc_single_simulationstate(c,v,eos);});
  }
// ***End*** Added by Emma
#endif
}

hdsim1D::hdsim1D
(const PhysicalGeometry1D& pg,
 const SimulationState1D& ss,
 const EquationOfState& eos,
 const VertexMotion& vm,
 const SourceTerm1D& force,
 const TimeStepFunction1D& tsf,
 const FluxCalculator1D& fc,
 const ExtensiveUpdater1D& eu,
 const CellUpdater1D& cu):
  pg_(pg),
  ss_(ss),
  eos_(eos),
  extensives_(calc_extensives(pg,ss,eos)),
  vm_(vm),
  force_(force),
  tsf_(tsf),
  fc_(fc),
  eu_(eu),
  cu_(cu),
  time_(0),
  cycle_(0),
  tracers_intensive_(vector<vector<double> >()),
  tracers_extensive_(vector<vector<double> >())
{
  spdlog::debug("hdsim1D initialisation completed");
}

void hdsim1D::TimeAdvance(void)
{
  spdlog::debug("begin time advance 1o iteration {0}, virtual time {1}",
		cycle_, time_);

  const vector<double> _VertexVelocity =
    serial_generate<size_t, double>
    (create_range<size_t>(0, ss_.getVertices().size()),
     [&](size_t i){return vm_(i,ss_.getVertices(), ss_.getCells());});

  spdlog::debug("Vertex velocity calculated");

  const double dt = tsf_(ss_,eos_);

  spdlog::debug("Time step calculated");

  const vector<Extensive> fluxes =
    fc_(ss_, _VertexVelocity, eos_, dt);

  spdlog::debug("Fluxes calculated");

  eu_
    (fluxes,
     pg_,
     ss_,
     dt,
     extensives_);

  spdlog::debug("Extensives updated");

  transform(extensives_.begin(),
	    extensives_.end(),
	    create_range<size_t>(0, extensives_.size()).begin(),
	    extensives_.begin(),
	    [&](const Extensive& e, size_t i)
	    {return e+dt*force_(ss_,i,fluxes,pg_,time_,dt);});

  spdlog::debug("source term calculated");

  transform(ss_.vertices_.begin(),
	    ss_.vertices_.end(),
	    _VertexVelocity.begin(),
	    ss_.vertices_.begin(),
	    [&](double x, double v)
	    {return x+v*dt;});

  spdlog::debug("Vertices updated");

  ss_.updateCells(cu_(pg_,
		      extensives_,
		      ss_,
		      eos_));

  spdlog::debug("cells updated");

  time_ += dt;
  cycle_++;

  spdlog::debug("time advance interation finished");
}

void hdsim1D::TimeAdvance2(void)
{
  const vector<double> mid_vertex_velocities =
    serial_generate<size_t, double>
    (create_range<size_t>(0, ss_.getVertices().size()),
     [&](size_t i){return vm_(i, ss_.getVertices(), ss_.getCells());});

  const double dt = tsf_(ss_, eos_);

  const vector<Extensive> mid_fluxes =
    fc_(ss_, mid_vertex_velocities, eos_, dt);

  vector<Extensive> mid_extensive = extensives_;

  eu_(mid_fluxes,
      pg_,
      ss_,
      dt/2,
      mid_extensive);

  transform(extensives_.begin(),
	    extensives_.end(),
	    create_range<size_t>(0, extensives_.size()).begin(),
	    extensives_.begin(),
	    [&](const Extensive& e, size_t i)
	    {return e+0.5*dt*force_(ss_,i,mid_extensive,pg_,time_,dt/2);});

  /*
  force_contribution(ss_,
		     force_, 
		     mid_fluxes,
		     pg_,
		     time_, 
		     dt/2,
		     mid_extensive);
  */

  SimulationState1D mid_state = ss_;
  transform(ss_.vertices_.begin(),
	    ss_.vertices_.end(),
	    mid_vertex_velocities.begin(),
	    ss_.vertices_.begin(),
	    [&](double x, double v)
	    {return x+dt*v;});
  /*
  mid_state.updateVertices
    (calc_new_vertices(mid_vertex_velocities,
		       dt,
		       ss_.getVertices()));
  */

  mid_state.updateCells
    (cu_(pg_,
	 mid_extensive,
	 mid_state,
	 eos_));

  const vector<double> _VertexVelocity = serial_generate<size_t, double>
    (create_range<size_t>(0, ss_.getVertices().size()),
     [&](size_t i){return vm_(i, ss_.getVertices(), ss_.getCells());});

  const vector<Extensive> fluxes =
    fc_(mid_state, _VertexVelocity, eos_, dt);

  eu_(fluxes,
      pg_,
      ss_,
      dt,
      extensives_);

  transform(extensives_.begin(),
	    extensives_.end(),
	    create_range<size_t>(0, extensives_.size()).begin(),
	    extensives_.begin(),
	    [&](const Extensive& e, size_t i)
	    {return e+dt*force_(mid_state, i, fluxes, pg_, time_, dt);});

  transform(ss_.vertices_.begin(),
	    ss_.vertices_.end(),
	    _VertexVelocity.begin(),
	    ss_.vertices_.begin(),
	    [&](double x, double v)
	    {return x+v*dt;});

  ss_.updateCells(cu_(pg_,
		      extensives_,
		      ss_,
		      eos_));

  time_ += dt;
  ++cycle_;
}

void hdsim1D::recalculateExtensives(void)
{
  extensives_ = calc_extensives
    (pg_,
     ss_,
     eos_);
}


// *** Start *** Added by Emma
void hdsim1D::recalculatePrimitives(void)
{
    ss_.updateCells(cu_(pg_,
		      extensives_,
		      ss_,
		      eos_));

}
// *** End *** Added by Emma

