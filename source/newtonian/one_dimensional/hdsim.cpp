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

namespace {
  vector<double> CalcVertexVelocities
  (const SimulationState1D& state, 
   VertexMotion const& vm)
  {
    vector<double> res(state.getVertices().size());
    for(size_t i = 0; i<state.getVertices().size();i++)
      res[i] = vm(i, state.getVertices(), state.getCells());

    return res;
  }

  vector<double> calc_new_vertices
  (const vector<double>& vv_,
   double dt,
   const vector<double>& vertices)
  {
    vector<double> res = vertices;
    for(size_t i=0;i<res.size();++i)
      res.at(i) += dt*vv_.at(i);
    return res;
  }
}

namespace {

  /*
  Extensive conserved2extensive
  (const Conserved& c)
  {
    Extensive res;
    res.mass = c.Mass;
    res.momentum = c.Momentum;
    res.energy = c.Energy;
    return res;
  }
  */
  
  void force_contribution
  (const SimulationState1D& state,
   const SourceTerm1D& force,
   const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   double t,
   double dt,
   vector<Extensive>& extensive)
  {
    for(size_t i=0;i<extensive.size();++i)
      extensive[i] +=
	dt*force(state, i, fluxes, pg, t, dt); 
  }
}

void hdsim1D::TimeAdvance(void)
{
  spdlog::debug("begin time advance 1o iteration {0}, virtual time {1}",
		cycle_, time_);

  /*
  const vector<double> _VertexVelocity = CalcVertexVelocities
    (ss_, vm_);
  */
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

  force_contribution(ss_,
		     force_, 
		     fluxes,
		     pg_,
		     time_, 
		     dt, 
		     extensives_);

  spdlog::debug("source term calculated");

  ss_.updateVertices(calc_new_vertices(_VertexVelocity,
				       dt,
				       ss_.getVertices()));

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
    CalcVertexVelocities(ss_, vm_);

  const double dt = tsf_(ss_, eos_);

  const vector<Extensive> mid_fluxes =
    fc_(ss_, mid_vertex_velocities, eos_, dt);

  vector<Extensive> mid_extensive = extensives_;

  eu_(mid_fluxes,
      pg_,
      ss_,
      dt/2,
      mid_extensive);
  
  force_contribution(ss_,
		     force_, 
		     mid_fluxes,
		     pg_,
		     time_, 
		     dt/2,
		     mid_extensive);

  SimulationState1D mid_state = ss_;
  mid_state.updateVertices
    (calc_new_vertices(mid_vertex_velocities,
		       dt,
		       ss_.getVertices()));

  /*
  const vector<Conserved> mid_intesive = 
    UpdateConservedIntensive
    (extensive2conserved(mid_extensive),
     mid_state.getVertices(),
     pg_);
  */
  mid_state.updateCells
    (cu_(pg_,
	 mid_extensive,
	 mid_state,
	 eos_));

  const vector<double> _VertexVelocity = CalcVertexVelocities
    (mid_state, vm_);

  const vector<Extensive> fluxes =
    fc_(mid_state, _VertexVelocity, eos_, dt);

  eu_(fluxes,
      pg_,
      ss_,
      dt,
      extensives_);

  force_contribution
    (mid_state,
     force_, 
     fluxes,
     pg_,
     time_, 
     dt,
     extensives_);

  ss_.updateVertices(calc_new_vertices(_VertexVelocity,
				       dt,
				       ss_.getVertices()));

  ss_.updateCells(cu_(pg_,
		      extensives_,
		      ss_,
		      eos_));

  time_ += dt;
  ++cycle_;
}
