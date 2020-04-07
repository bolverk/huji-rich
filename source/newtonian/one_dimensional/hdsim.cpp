#include <cmath>
#include <algorithm>
#include "hdsim.hpp"
#include "../../misc/universal_error.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"
#include "spdlog/spdlog.h"

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

// External functions

namespace {

  double GetVolume
  (const vector<double>& v, 
   const PhysicalGeometry1D& pg,
   size_t i)
  {
    return pg.calcVolume(v.at(i+1)) 
      - pg.calcVolume(v.at(i));
  }

  vector<Extensive> calc_extensives
  (const PhysicalGeometry1D& pg,
   const SimulationState1D& ss,
   const EquationOfState& eos)
  {
    const vector<double>& vertices = ss.getVertices();
    const vector<ComputationalCell>& cells = ss.getCells();
    vector<Extensive> res(ss.getCells().size());
    for(size_t i=0;i<res.size();++i){
      const double volume = GetVolume(vertices,
				      pg,
				      i);
      res.at(i).mass = cells.at(i).density*volume;
      res.at(i).momentum = res.at(i).mass*cells.at(i).velocity;
      const double kinetic_specific_energy =
	0.5*pow(abs(cells.at(i).velocity),2);
      const double thermal_specific_energy =
	eos.dp2e(cells.at(i).density, cells.at(i).pressure);
      res.at(i).energy = res.at(i).mass*
	(kinetic_specific_energy+thermal_specific_energy);
      for(size_t j=0;j<cells.at(0).tracers.size();++j){
	res.at(i).tracers.push_back
	  (cells.at(i).tracers.at(j)*res.at(i).mass);
      }
    }
    
    return res;
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
   double t,
   double dt,
   vector<Extensive>& extensive)
  {
    for(size_t i=0;i<extensive.size();++i)
      extensive[i] +=
	dt*force(state, i, t, dt); 
  }
}

void hdsim1D::TimeAdvance(void)
{
  spdlog::debug("begin time advance 1o iteration {0}, virtual time {1}",
		cycle_, time_);
  
  const vector<double> _VertexVelocity = CalcVertexVelocities
    (ss_, vm_);

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
		     force_, time_, dt, 
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
		     force_, time_, dt/2,
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
     force_, time_, dt,
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
