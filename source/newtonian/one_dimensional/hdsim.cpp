#include <cmath>
#include <algorithm>
#include "hdsim.hpp"
#include "../../misc/universal_error.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/utils.hpp"
#include "../../misc/lazy_list.hpp"

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

namespace{

  vector<Conserved> extensive2conserved
    (const vector<Extensive>& extensives)
  {
    vector<Conserved> res(extensives.size());
    for(size_t i=0;i<res.size();++i){
      res.at(i).Mass = extensives.at(i).mass;
      res.at(i).Momentum = extensives.at(i).momentum;
      res.at(i).Energy = extensives.at(i).energy;
    }
    return res;
  }

  vector<Extensive> conserved2extensive
    (const vector<Conserved>& conserved)
  {
    vector<Extensive> res(conserved.size());
    for(size_t i=0;i<res.size();++i){
      res.at(i).mass = conserved.at(i).Mass;
      res.at(i).momentum = conserved.at(i).Momentum;
      res.at(i).energy = conserved.at(i).Energy;
    }
    return res;
  }
  
  vector<Primitive> cc2primitives
    (const vector<ComputationalCell>& ccs,
     const EquationOfState& eos)
  {
    vector<Primitive> res(ccs.size());
    for(size_t i=0;i<res.size();++i){
      res.at(i).Density = ccs.at(i).density;
      res.at(i).Pressure = ccs.at(i).pressure;
      res.at(i).Velocity = ccs.at(i).velocity;
      res.at(i).Energy = eos.dp2e
	(ccs.at(i).density, ccs.at(i).pressure);
      res.at(i).SoundSpeed = eos.dp2c
	(ccs.at(i).density, ccs.at(i).pressure);
    }
    return res;
  }

  /*
  vector<ComputationalCell> primitives2cc
    (const vector<Primitive>& primitives)
  {
    vector<ComputationalCell> res(primitives.size());
    for(size_t i=0;i<res.size();++i){
      res.at(i).density = primitives.at(i).Density;
      res.at(i).pressure = primitives.at(i).Pressure;
      res.at(i).velocity = primitives.at(i).Velocity;
    }
    return res;
  }
  */
}

const SimulationState1D& hdsim1D::getState(void) const
{
  return ss_;
}

/*
const vector<Primitive> hdsim1D::getCells(void) const
{
  return cc2primitives(ss_.getCells(),eos_);
  }

void hdsim1D::setCells(const vector<Primitive>& primitives)
{
  //  _Cells = primitives;
  ss_.updateCells(primitives2cc(primitives));
}
*/

//double hdsim1D::GetVertexPosition(size_t i) const
//{
  //  return _Vertices[i];
//  return ss_.getVertices().at(i);
//}

// External functions

namespace {
  vector<Conserved> CalcConservedIntensive(const vector<Primitive>& p)
  {
    vector<Conserved> res(p.size());
    for(size_t i=0;i<p.size();i++){
      res[i] = Primitive2Conserved(p[i]);
    }
    return res;
  }

  double GetVolume
  (const vector<double>& v, 
   const PhysicalGeometry1D& pg,
   size_t i)
  {
    return pg.calcVolume(v.at(i+1)) 
      - pg.calcVolume(v.at(i));
  }

  vector<Conserved> CalcConservedExtensive
  (const PhysicalGeometry1D& pg,
   const vector<Conserved>& ci, 
   const vector<double>& v)
  {
    vector<Conserved> res(ci.size());
    for(size_t i=0;i<ci.size();i++){
      res[i] = GetVolume(v, pg, i)*ci[i];
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
  extensives_
  (conserved2extensive
   (CalcConservedExtensive
    (pg_,
     CalcConservedIntensive(cc2primitives(ss_.getCells(),eos)),
     ss_.getVertices()))),
  vm_(vm),
  force_(force),
  tsf_(tsf),
  fc_(fc),
  eu_(eu),
  cu_(cu),
  time_(0),
  cycle_(0),
  tracers_intensive_(vector<vector<double> >()),
  tracers_extensive_(vector<vector<double> >()) {}

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
  (const vector<double> vv_,
   double dt,
   const vector<double>& vertices)
  {
    vector<double> res = vertices;
    for(size_t i=0;i<res.size();++i)
      res.at(i) += dt*vv_.at(i);
    return res;
  }

  vector<Conserved> UpdateConservedIntensive
  (const vector<Conserved>& ConservedExtensive, 
   const vector<double>& Vertices,
   const PhysicalGeometry1D& pg)
  {
    vector<Conserved> res(ConservedExtensive.size());
    for(size_t i=0;i<ConservedExtensive.size();i++){
      res[i] = ConservedExtensive[i] / 
	GetVolume(Vertices, pg, i);
    }
    return res;
  }
}

namespace {

  Extensive conserved2extensive
  (const Conserved& c)
  {
    Extensive res;
    res.mass = c.Mass;
    res.momentum = c.Momentum;
    res.energy = c.Energy;
    return res;
  }
  
  void force_contribution
  (const SimulationState1D& state,
   const SourceTerm1D& force,
   double t,
   double dt,
   vector<Extensive>& extensive)
  {
    for(size_t i=0;i<extensive.size();++i)
      extensive[i] +=
	dt*conserved2extensive(force(state,
				     i, t, dt)); 
  }
}

void hdsim1D::TimeAdvance(void)
{
  const vector<double> _VertexVelocity = CalcVertexVelocities
    (ss_, vm_);

  //  const double dt = _cfl*MaxTimeStep(ss_.getVertices(), getCells());
  const double dt = tsf_(ss_,eos_);

  const vector<Extensive> fluxes =
    fc_(ss_, _VertexVelocity, eos_, dt);

  eu_
    (fluxes,
     pg_,
     ss_,
     dt,
     extensives_);

  force_contribution(ss_,
		     force_, time_, dt, 
		     extensives_);

  //  MoveVertices(_VertexVelocity, dt, ss_.getVertices());
  ss_.updateVertices(calc_new_vertices(_VertexVelocity,
				       dt,
				       ss_.getVertices()));

  const vector<Conserved> _ConservedIntensive = UpdateConservedIntensive
    (extensive2conserved(extensives_),
     ss_.getVertices(), pg_);

  ss_.updateCells(cu_(_ConservedIntensive,
		      extensive2conserved(extensives_),
		      ss_.getCells(),
		      eos_));

  time_ += dt;
  cycle_++;
}

void hdsim1D::TimeAdvance2(void)
{
  const vector<double> mid_vertex_velocities = 
    CalcVertexVelocities(ss_, vm_);

  //  const double dt = _cfl*MaxTimeStep(ss_.getVertices(), getCells());
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

  const vector<Conserved> mid_intesive = 
    UpdateConservedIntensive
    (extensive2conserved(mid_extensive),
     mid_state.getVertices(),
     pg_);
  mid_state.updateCells
    (cu_(mid_intesive,
	 extensive2conserved(mid_extensive),
	 mid_state.getCells(),
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

  const vector<Conserved> _ConservedIntensive = UpdateConservedIntensive
    (extensive2conserved(extensives_),
     ss_.getVertices(), pg_);

  ss_.updateCells(cu_(_ConservedIntensive,
				    extensive2conserved(extensives_),
				    ss_.getCells(),
				    eos_));

  time_ += dt;
  ++cycle_;
}
