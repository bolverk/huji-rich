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

double hdsim1D::GetCellCenter(size_t index) const
{
  return 0.5*(ss_.getVertices().at(index)+
	      ss_.getVertices().at(index+1));
}

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
}

const vector<Primitive> hdsim1D::getCells(void) const
{
  //  return _Cells;
  return cc2primitives(ss_.getCells(),_eos);
}

void hdsim1D::setCells(const vector<Primitive>& primitives)
{
  //  _Cells = primitives;
  ss_.updateCells(primitives2cc(primitives));
}

int hdsim1D::GetVertexNo(void) const
{
  //  return static_cast<int>(_Vertices.size());
  return static_cast<int>(ss_.getVertices().size());
}

double hdsim1D::GetVertexPosition(size_t i) const
{
  //  return _Vertices[i];
  return ss_.getVertices().at(i);
}

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
 const vector<double>& vertices,
 const SpatialReconstruction1D& Interpolation,
 const SpatialDistribution1D& density,
 const SpatialDistribution1D& pressure,
 const SpatialDistribution1D& paravelocity,
 const SpatialDistribution1D& perpvelocity,
 const EquationOfState& eos,
 const RiemannSolver& rs,
 const VertexMotion& vm,
 const BoundaryConditions1D& bc,
 const SourceTerm1D& force,
 const TimeStepFunction1D& tsf,
 const FluxCalculator1D& fc,
 const ExtensiveUpdater1D& eu,
 const CellUpdater1D& cu):
  pg_(pg),
  ss_(vertices,
      density,
      pressure,
      paravelocity,
      perpvelocity,
      vector<pair<string, const SpatialDistribution1D*> >(),
      vector<pair<string, const BoolSpatialDistribution* > >()),
  _eos(eos), 
  _ConservedExtensive
  (conserved2extensive
   (CalcConservedExtensive
    (pg_,
     CalcConservedIntensive(cc2primitives(ss_.getCells(),eos)),
     ss_.getVertices()))),
  _Interpolation(Interpolation),
  _rs(rs), _vm(vm), _bc(bc), 
  force_(force),
  tsf_(tsf),
  fc_(fc),
  eu_(eu),
  cu_(cu),
  time_(0),
  cycle_(0),
  tracers_intensive_(vector<vector<double> >()),
  tracers_extensive_(vector<vector<double> >()),
  cold_flows_() {}

namespace {
  vector<double> CalcVertexVelocities
  (vector<double> const& Vertices, 
   vector<Primitive> const& Cells,
   VertexMotion const& vm)
  {
    vector<double> res(Vertices.size());
    for(size_t i = 0; i<Vertices.size();i++)
      res[i] = vm.CalcVelocity(int(i), Vertices, Cells);

    return res;
  }

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
    __attribute__((noreturn))
#endif
  void riemann_solver_rethrow
  (Primitive const& left,
   Primitive const& right,
   size_t idx, 
   double pos,
   double vertex_velocity,
   UniversalError& eo)
  {
    eo.AddEntry("riemann solver stage data starts here",0);
    eo.AddEntry("left density",left.Density);
    eo.AddEntry("left pressure",left.Pressure);
    eo.AddEntry("left x velocity",left.Velocity.x);
    eo.AddEntry("left y velocity",left.Velocity.y);
    eo.AddEntry("left sound speed",left.SoundSpeed);
    eo.AddEntry("right density",right.Density);
    eo.AddEntry("right pressure",right.Pressure);
    eo.AddEntry("right x velocity",right.Velocity.x);
    eo.AddEntry("right y velocity",right.Velocity.y);
    eo.AddEntry("right sound speed",right.SoundSpeed);
    eo.AddEntry("right energy",right.Energy);
    eo.AddEntry("interface index",static_cast<double>(idx));
    eo.AddEntry("interface position",pos);
    eo.AddEntry("interface velocity",vertex_velocity);
    throw eo;
  }

  vector<Conserved> SolveRiemannProblems
  (vector<double> const& Vertices, 
   vector<Primitive> const& Cells,
   SpatialReconstruction1D const& Interpolation,
   vector<double> const& VertexVelocity,
   RiemannSolver const& rs,
   BoundaryConditions1D const& bc,
   double dt)
  {
    vector<Conserved> res(Vertices.size());
    for(size_t i = 1;i<Vertices.size()-1; i++){
      const Primitive left = Interpolation.InterpState
	(Vertices, Cells, VertexVelocity[i], i, 0,dt);
      const Primitive right = Interpolation.InterpState
	(Vertices, Cells, VertexVelocity[i], i, 1,dt);
      try{
	res[i] = rs(left, right, VertexVelocity[i]);
      }
      catch(UniversalError& eo){
	riemann_solver_rethrow(left,
			       right,
			       i, 
			       Vertices[i],
			       VertexVelocity[i],
			       eo);
      }
    }
    res[0] = bc.CalcFlux(Vertices, Cells, rs, 
			 VertexVelocity,0);
    res[Vertices.size()-1] = 
      bc.CalcFlux(Vertices, Cells, rs, 
		  VertexVelocity, static_cast<int>(Vertices.size())-1);
    return res;
  }

  void MoveVertices(vector<double> const& VertexVelocity,
		    double dt, vector<double>& Vertices)
  {
    for(size_t i=0;i<Vertices.size();i++){
      Vertices[i] += dt*VertexVelocity[i];
    }
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
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   const SourceTerm1D& force,
   double t,
   double dt,
   vector<Extensive>& extensive)
  {
    for(size_t i=0;i<extensive.size();++i)
      extensive[i] +=
	dt*conserved2extensive(force(vertices,
				     cells,
				     i, t, dt)); 
  }
}

void hdsim1D::TimeAdvance(void)
{
  cold_flows_.initializeEntropies(ss_.getVertices(),
				  getCells(),
				  _eos);

  const vector<double> _VertexVelocity = CalcVertexVelocities
    (ss_.getVertices(), getCells(), _vm);

  //  const double dt = _cfl*MaxTimeStep(ss_.getVertices(), getCells());
  const double dt = tsf_(ss_,_eos);

  /*
  const vector<Conserved> _Fluxes = SolveRiemannProblems
    (ss_.getVertices(), getCells(), _Interpolation, _VertexVelocity,
     _rs, _bc, dt);
  */
  const vector<Extensive> fluxes =
    fc_(ss_, _VertexVelocity, _eos, dt);

  cold_flows_.advanceEntropies
    (extensive2conserved(fluxes),
     extensive2conserved(_ConservedExtensive),
     dt);

  eu_
    (fluxes,
     pg_,
     ss_,
     dt,
     _ConservedExtensive);

  force_contribution(ss_.getVertices(), getCells(),
		     force_, time_, dt, 
		     _ConservedExtensive);

  //  MoveVertices(_VertexVelocity, dt, ss_.getVertices());
  ss_.updateVertices(calc_new_vertices(_VertexVelocity,
				       dt,
				       ss_.getVertices()));

  const vector<Conserved> _ConservedIntensive = UpdateConservedIntensive
    (extensive2conserved(_ConservedExtensive),
     ss_.getVertices(), pg_);

  if(cold_flows_.is_active())
    setCells(cold_flows_.retrieveAllPrimitive
	     (_ConservedIntensive,
	      extensive2conserved(_ConservedExtensive),
	      _eos));
  else
    setCells(cu_(_ConservedIntensive,
		 extensive2conserved(_ConservedExtensive),
		 getCells(),
		 _eos));

  time_ += dt;
  cycle_++;
}

void hdsim1D::TimeAdvance2(void)
{
  cold_flows_.initializeEntropies(ss_.getVertices(),
				  getCells(),
				  _eos);

  const vector<double> mid_vertex_velocities = 
    CalcVertexVelocities(ss_.getVertices(), getCells(), _vm);

  //  const double dt = _cfl*MaxTimeStep(ss_.getVertices(), getCells());
  const double dt = tsf_(ss_, _eos);

  const vector<Conserved> mid_fluxes = 
    SolveRiemannProblems(ss_.getVertices(), getCells(), _Interpolation, 
			 mid_vertex_velocities,
			 _rs, _bc, dt);

  cold_flows_.advanceEntropies
    (mid_fluxes, extensive2conserved(_ConservedExtensive), dt/2);

  vector<Extensive> mid_extensive = _ConservedExtensive;

  eu_(conserved2extensive(mid_fluxes),
      pg_,
      ss_,
      dt/2,
      mid_extensive);
  
  force_contribution(ss_.getVertices(), getCells(),
		     force_, time_, dt/2,
		     mid_extensive);
  vector<double> mid_vertices = ss_.getVertices();
  MoveVertices(mid_vertex_velocities, 
	       dt/2, mid_vertices);

  const vector<Conserved> mid_intesive = 
    UpdateConservedIntensive
    (extensive2conserved(mid_extensive),
     mid_vertices,
     pg_);

  const vector<Primitive> mid_cells = 
    cold_flows_.retrieveAllPrimitive(mid_intesive,
				     extensive2conserved(mid_extensive),
				     _eos);

  cold_flows_.initializeEntropies(ss_.getVertices(),
				  getCells(),
				  _eos);

  const vector<double> _VertexVelocity = CalcVertexVelocities
    (mid_vertices, mid_cells, _vm);

  const vector<Conserved> _Fluxes = SolveRiemannProblems
    (mid_vertices, mid_cells, _Interpolation, _VertexVelocity,
     _rs, _bc, dt);

  cold_flows_.advanceEntropies
    (_Fluxes,
     CalcConservedIntensive(getCells()),
     dt);

  eu_(conserved2extensive(_Fluxes),
      pg_,
      ss_,
      dt,
      _ConservedExtensive);

  force_contribution(mid_vertices, mid_cells,
		     force_, time_, dt,
		     _ConservedExtensive);

  //  MoveVertices(_VertexVelocity, dt, _Vertices);
  ss_.updateVertices(calc_new_vertices(_VertexVelocity,
				       dt,
				       ss_.getVertices()));

  const vector<Conserved> _ConservedIntensive = UpdateConservedIntensive
    (extensive2conserved(_ConservedExtensive),
     ss_.getVertices(), pg_);

  setCells(cold_flows_.retrieveAllPrimitive
	   (_ConservedIntensive,
	    extensive2conserved(_ConservedExtensive),
	    _eos));

  time_ += dt;
  ++cycle_;
}

ColdFlows::ColdFlows(void):
  active_(false),
  threshold_(0),
  entropies_() {}

void ColdFlows::activate(double threshold)
{
  active_ = true;
  assert(threshold>0);
  threshold_ = threshold;
}

namespace {
  class EntropyCalculator: public LazyList<double>
  {
  public:

    EntropyCalculator(const vector<double>& grid,
		      const vector<Primitive>& cells,
		      const EquationOfState& eos):
      grid_(grid), cells_(cells), eos_(eos) {}

    size_t size(void) const
    {
      return cells_.size();
    }

    double operator[](size_t i) const
    {
      const double volume = grid_[i+1] - grid_[i];
      const double mass = volume*cells_[i].Density;
      return mass*eos_.dp2s(cells_[i].Density,
			    cells_[i].Pressure);
    }

  private:
    const vector<double>& grid_;
    const vector<Primitive>& cells_;
    const EquationOfState& eos_;
  };
}

void ColdFlows::initializeEntropies(const vector<double>& grid,
					     const vector<Primitive>& cells,
					     const EquationOfState& eos)
{
  if(active_)
    entropies_ = serial_generate(EntropyCalculator(grid,cells,eos));
}

void ColdFlows::advanceEntropies(const vector<Conserved>& fluxes,
					  const vector<Conserved>& extensive,
					  double dt)
{
  if(active_){
    vector<double> specific_entropies(extensive.size());
    for(size_t i=0;i<extensive.size();++i)
      specific_entropies[i] = entropies_[i]/extensive[i].Mass;
    vector<double> entropy_fluxes(fluxes.size());
    entropy_fluxes[0] = specific_entropies[0]*fluxes[0].Mass;
    for(size_t i=1;i<extensive.size();++i)
      entropy_fluxes[i] = fluxes[i].Mass*specific_entropies[fluxes[i].Mass>0 ? i-1 : i];
    entropy_fluxes[entropy_fluxes.size()-1] = specific_entropies.back()*fluxes.back().Mass;
    vector<double> entropy_difference(extensive.size());
    for(size_t i=0;i<extensive.size();++i)
      entropy_difference[i] = dt*(entropy_fluxes[i] - entropy_fluxes[i+1]);
    vector<double> old_entropies = entropies_;
    for(size_t i=0;i<entropies_.size();++i)
      entropies_[i] = old_entropies[i] + entropy_difference[i];

    for(size_t i=0;i<entropies_.size();++i)
      assert(entropies_[i]>0 && "Entropy is negative");
  }
}

namespace {
  Primitive retrieve_single_primitive(const Conserved& intensive,
				      const EquationOfState& eos,
				      double specific_entropy,
				      double thres,
				      bool active)
  {
    const double density = intensive.Mass;
    const Vector2D velocity = intensive.Momentum/intensive.Mass;
    const double kinetic_energy = 0.5*(pow(abs(velocity),2));
    const double thermal_energy = intensive.Energy/intensive.Mass -
      kinetic_energy;
    const double pressure = (active&&(thermal_energy/kinetic_energy<thres)) ?
      eos.sd2p(specific_entropy, density) :
      eos.de2p(density, thermal_energy);
    const double sound_speed = eos.dp2c(density, pressure);
    return Primitive(density, pressure, velocity,
		     thermal_energy, sound_speed);
  }

  class PrimitiveRetriever: public LazyList<Primitive>
  {
  public:

    PrimitiveRetriever(const vector<Conserved>& intensive,
		       const vector<Conserved>& extensive,
		       const EquationOfState& eos,
		       const vector<double>& entropies,
		       double threshold,
		       bool active):
      intensive_(intensive),
      extensive_(extensive),
      eos_(eos),
      entropies_(entropies),
      threshold_(threshold),
      active_(active) {}

    size_t size(void) const
    {
      return intensive_.size();
    }

    Primitive operator[](size_t i) const
    {
      return retrieve_single_primitive
	(intensive_.at(i),eos_,
	 (active_ ? entropies_.at(i)/extensive_.at(i).Mass : 0),
	 threshold_, active_);
    }
    
  private:
    const vector<Conserved>& intensive_;
    const vector<Conserved>& extensive_;
    const EquationOfState& eos_;
    const vector<double> entropies_;
    const double threshold_;
    const bool active_;
  };
}

vector<Primitive> ColdFlows::retrieveAllPrimitive
(const vector<Conserved>& intensive,
 const vector<Conserved>& extensive,
 const EquationOfState& eos) const
{
  return serial_generate(PrimitiveRetriever(intensive,
					    extensive,
					    eos,
					    entropies_,
					    threshold_,
					    active_));
}

bool ColdFlows::is_active(void) const
{
    return active_;
}

void hdsim1D::enableColdFlows(double thres)
{
  cold_flows_.activate(thres);
}
