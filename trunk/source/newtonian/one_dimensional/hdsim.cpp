#include <cmath>
#include <algorithm>
#include "hdsim.hpp"
#include "../../misc/universal_error.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../common/hydrodynamics.hpp"
#include "../../misc/utils.hpp"

using namespace std;

HydroSnapshot1D::HydroSnapshot1D
(vector<double> const& redges,
 vector<Primitive> const& rcells,
 vector<Conserved> const& rintensive,
 vector<Conserved> const& rextensive):
  edges(redges), cells(rcells), 
  intensive(rintensive),
  extensive(rextensive) {}

namespace {
  vector<Primitive> InitialiseCells
  (vector<double> const& vertices,
   SpatialDistribution1D const& density,
   SpatialDistribution1D const& pressure,
   SpatialDistribution1D const& paravelocity,
   SpatialDistribution1D const& perpvelocity,
   EquationOfState const& eos)
  {
    vector<Primitive> res(vertices.size()-1);
    for(size_t i = 0; i<vertices.size() - 1; i++){
      const double r = 0.5*(vertices[i] + vertices[i+1]);
      const double d = density.EvalAt(r);
      const double p = pressure.EvalAt(r);
      const Vector2D v(paravelocity.EvalAt(r),
		       perpvelocity.EvalAt(r));
      res[i] = CalcPrimitive(d, p, v, eos);
    }
    return res;
  }
}

// Diagnostics

double hdsim1D::GetCellCenter(size_t index) const
{
  return 0.5*(_Vertices[index]+
	      _Vertices[index+1]);
}

double hdsim1D::GetTime(void) const
{
  return time_;
}

int hdsim1D::GetCycle(void) const
{
  return cycle_;
}

vector<Conserved> const& hdsim1D::getFluxes(void) const
{
  return _Fluxes;
}

int hdsim1D::GetVertexNo(void) const
{
  return (int)_Vertices.size();
}

double hdsim1D::GetVertexPosition(size_t i) const
{
  return _Vertices[i];
}

int hdsim1D::GetCellNo(void) const
{
  return (int)_Cells.size();
}

Primitive hdsim1D::GetCell(size_t i) const
{
  return _Cells[i];
}

// External functions

namespace {
  vector<Conserved> CalcConservedIntensive(vector<Primitive> p)
  {
    vector<Conserved> res(p.size());
    for(size_t i=0;i<p.size();i++){
      res[i] = Primitive2Conserved(p[i]);
    }
    return res;
  }

  double GetVolume(vector<double> v, size_t i)
  {
    return v[i+1] - v[i];
  }

  vector<Conserved> CalcConservedExtensive
  (vector<Conserved> const& ci, vector<double> const& v)
  {
    vector<Conserved> res(ci.size());
    for(size_t i=0;i<ci.size();i++){
      res[i] = GetVolume(v, i)*ci[i];
    }
    return res;
  }
}

hdsim1D::hdsim1D
(vector<double> const& vertices,
 SpatialReconstruction1D const& Interpolation,
 SpatialDistribution1D const& density,
 SpatialDistribution1D const& pressure,
 SpatialDistribution1D const& paravelocity,
 SpatialDistribution1D const& perpvelocity,
 EquationOfState const& eos,
 RiemannSolver const& rs,
 VertexMotion const& vm,
 BoundaryConditions1D const& bc,
 ExternalForces1D const& force):
  _Vertices(vertices), _eos(eos), 
  _Cells(InitialiseCells(vertices, density, pressure,
			 paravelocity, perpvelocity, eos)),
  _Fluxes(vector<Conserved>(vertices.size())),
  _VertexVelocity(vector<double>()),
  _ConservedIntensive(CalcConservedIntensive(_Cells)),
  _ConservedExtensive(CalcConservedExtensive
		      (_ConservedIntensive, _Vertices)),
  _Interpolation(Interpolation),
  _rs(rs), _vm(vm), _bc(bc), 
  force_(force), _cfl(1./3.), time_(0), cycle_(0),
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

  __attribute__((noreturn)) void riemann_solver_rethrow
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
    eo.AddEntry("interface index",(double)idx);
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
	(Vertices, Cells, VertexVelocity[i], int(i), 0,dt);
      const Primitive right = Interpolation.InterpState
	(Vertices, Cells, VertexVelocity[i], int(i), 1,dt);
      try{
	res[i] = rs.Solve(left, right, VertexVelocity[i]);
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
		  VertexVelocity, (int)Vertices.size()-1);
    return res;
  }

  double MaxTimeStepForCell(double width, Primitive const& p)
  {
    return width/(p.SoundSpeed+abs(p.Velocity.x));
  }

  double MaxTimeStep(vector<double> const& Vertices,
		     vector<Primitive> const& Cells)
  {
    double res = MaxTimeStepForCell(Vertices[1]-Vertices[0], Cells[0]);
    for(size_t i=1;i<Vertices.size()-1;i++){
      res = min(res,MaxTimeStepForCell
		(Vertices[i+1]-Vertices[i],Cells[i]));
    }
    return res;
  }

  void UpdateConservedExtensive
  (vector<Conserved> const& Fluxes, double dt,
   vector<Conserved>& ConservedExtensive)
  {
    for(size_t i = 0; i<ConservedExtensive.size(); i++){
      ConservedExtensive[i] += dt*Fluxes[i];
      ConservedExtensive[i] -= dt*Fluxes[i+1];
    }
  }

  void MoveVertices(vector<double> const& VertexVelocity,
		    double dt, vector<double>& Vertices)
  {
    for(size_t i=0;i<Vertices.size();i++){
      Vertices[i] += dt*VertexVelocity[i];
    }
  }

  vector<Conserved> UpdateConservedIntensive
  (vector<Conserved> const& ConservedExtensive, 
   vector<double> const& Vertices)
  {
    vector<Conserved> res(ConservedExtensive.size());
    for(size_t i=0;i<ConservedExtensive.size();i++){
      res[i] = ConservedExtensive[i] / 
	GetVolume(Vertices, i);
    }
    return res;
  }

  vector<Primitive> UpdatePrimitives
  (vector<Conserved> const& ConservedIntensive,
   EquationOfState const& eos)
  {
    vector<Primitive> res(ConservedIntensive.size());
    for(size_t i=0;i<ConservedIntensive.size();i++)
      res[i] = Conserved2Primitive(ConservedIntensive[i], eos);
    return res;
  }
}

void hdsim1D::overrideCFL(double cfl)
{
  _cfl = cfl;
}

namespace {
  void force_contribution
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   ExternalForces1D const& force,
   double t,
   double dt,
   vector<Conserved>& extensive)
  {
    for(size_t i=0;i<extensive.size();++i)
      extensive[i] +=
	dt*force.calc(vertices,
		      cells,
		      int(i), t, dt); 
  }
}

void hdsim1D::TimeAdvance(void)
{
  cold_flows_.initializeEntropies(_Vertices,
				  _Cells,
				  _eos);

  _VertexVelocity = CalcVertexVelocities
    (_Vertices, _Cells, _vm);

  const double dt = _cfl*MaxTimeStep(_Vertices, _Cells);

  _Fluxes = SolveRiemannProblems
    (_Vertices, _Cells, _Interpolation, _VertexVelocity,
     _rs, _bc, dt);

  cold_flows_.advanceEntropies(_Fluxes,_ConservedExtensive,dt);

  UpdateConservedExtensive(_Fluxes, dt, _ConservedExtensive);

  force_contribution(_Vertices, _Cells,
		     force_, time_, dt, 
		     _ConservedExtensive);

  MoveVertices(_VertexVelocity, dt, _Vertices);

  _ConservedIntensive = UpdateConservedIntensive
    (_ConservedExtensive, _Vertices);

  //  _Cells = UpdatePrimitives(_ConservedIntensive, _eos);
  _Cells = cold_flows_.retrieveAllPrimitive(_ConservedIntensive,
					    _ConservedExtensive,
					    _eos);

  time_ += dt;
  cycle_++;
}

void hdsim1D::TimeAdvance2(void)
{
  cold_flows_.initializeEntropies(_Vertices,
				  _Cells,
				  _eos);

  const vector<double> mid_vertex_velocities = 
    CalcVertexVelocities(_Vertices, _Cells, _vm);

  const double dt = _cfl*MaxTimeStep(_Vertices, _Cells);

  const vector<Conserved> mid_fluxes = 
    SolveRiemannProblems(_Vertices, _Cells, _Interpolation, 
			 mid_vertex_velocities,
			 _rs, _bc, dt);

  cold_flows_.advanceEntropies(mid_fluxes, _ConservedExtensive, dt/2);

  vector<Conserved> mid_extensive = _ConservedExtensive;
  UpdateConservedExtensive(mid_fluxes,dt/2,mid_extensive);
  force_contribution(_Vertices, _Cells,
		     force_, time_, dt/2,
		     mid_extensive);
  vector<double> mid_vertices = _Vertices;
  MoveVertices(mid_vertex_velocities, 
	       dt/2, mid_vertices);

  const vector<Conserved> mid_intesive = 
    UpdateConservedIntensive(mid_extensive, mid_vertices);

  const vector<Primitive> mid_cells = 
    cold_flows_.retrieveAllPrimitive(mid_intesive,
				     mid_extensive,
				     _eos);

  cold_flows_.initializeEntropies(_Vertices,
				  _Cells,
				  _eos);

  _VertexVelocity = CalcVertexVelocities
    (mid_vertices, mid_cells, _vm);

  _Fluxes = SolveRiemannProblems
    (mid_vertices, mid_cells, _Interpolation, _VertexVelocity,
     _rs, _bc, dt);

  cold_flows_.advanceEntropies(_Fluxes,_ConservedIntensive,dt);

  UpdateConservedExtensive(_Fluxes, dt, _ConservedExtensive);

  force_contribution(mid_vertices, mid_cells,
		     force_, time_, dt,
		     _ConservedExtensive);

  MoveVertices(_VertexVelocity, dt, _Vertices);

  _ConservedIntensive = UpdateConservedIntensive
    (_ConservedExtensive, _Vertices);

  _Cells = cold_flows_.retrieveAllPrimitive(_ConservedIntensive,
					    _ConservedExtensive,
					    _eos);

  time_ += dt;
  ++cycle_;
}

namespace{
  HydroSnapshot1D time_advance_1st_order
    (const HydroSnapshot1D& old,
     const VertexMotion& vm,
     const SpatialReconstruction1D& sr,
     const RiemannSolver& rs,
     const BoundaryConditions1D& bc,
     const EquationOfState& eos,
     const ExternalForces1D& force,
     double t, double dt)
  {
    const vector<double> edge_velocity = CalcVertexVelocities
      (old.edges,old.cells,vm);

    const vector<Conserved> fluxes = SolveRiemannProblems
      (old.edges,old.cells,sr,edge_velocity,rs,bc,dt);

    vector<Conserved> extensive = old.extensive;
    UpdateConservedExtensive(fluxes,dt,extensive);
    force_contribution(old.edges,
		       old.cells,
		       force,
		       t, dt,
		       extensive);
  
    vector<double> edges = old.edges;
    MoveVertices(edge_velocity,dt,edges);

    const vector<Conserved> intensive = UpdateConservedIntensive(extensive,edges);

    const vector<Primitive> cells = UpdatePrimitives(intensive,eos);

    return HydroSnapshot1D(edges,cells,intensive,extensive);
  }

  HydroSnapshot1D time_advance_2nd_order
    (HydroSnapshot1D const& old,
     VertexMotion const& vm,
     SpatialReconstruction1D const& sr,
     RiemannSolver const& rs,
     BoundaryConditions1D const& bc,
     EquationOfState const& eos,
     ExternalForces1D const& force,
     double t, double dt)
  {
    const HydroSnapshot1D mid = time_advance_1st_order
      (old,vm,sr,rs,bc,eos,force,t,dt/2);

    const vector<double> edge_velocity = CalcVertexVelocities
      (mid.edges,mid.cells,vm);

    const vector<Conserved> fluxes = SolveRiemannProblems
      (mid.edges,mid.cells,sr,edge_velocity,rs,bc,dt);

    vector<Conserved> new_extensive = old.extensive;
    UpdateConservedExtensive(fluxes,dt,new_extensive);
    force_contribution(mid.edges,
		       mid.cells,
		       force,
		       t+dt/2, dt,
		       new_extensive);

    vector<double> new_edges = old.edges;
    MoveVertices(edge_velocity,dt,new_edges);

    vector<Conserved> new_intensive = 
      UpdateConservedIntensive(new_extensive,new_edges);

    vector<Primitive> new_cells = 
      UpdatePrimitives(new_intensive,eos);

    return HydroSnapshot1D(new_edges,
			   new_cells,
			   new_intensive,
			   new_extensive);
  }
}

void hdsim1D::TimeAdvanceRK(int order)
{
  const double dt = _cfl*MaxTimeStep(_Vertices, _Cells);
  if(1==order){
    HydroSnapshot1D res = time_advance_1st_order
      (HydroSnapshot1D(_Vertices,_Cells,
		       _ConservedIntensive,
		       _ConservedExtensive),
       _vm, _Interpolation, _rs, _bc, _eos, 
       force_, time_, dt);
    _Vertices = res.edges;
    _Cells = res.cells;
    _ConservedIntensive = res.intensive;
    _ConservedExtensive = res.extensive;
  }
  else if(2==order){
    HydroSnapshot1D res = time_advance_2nd_order
      (HydroSnapshot1D(_Vertices,_Cells,
		       _ConservedIntensive,
		       _ConservedExtensive),
       _vm, _Interpolation, _rs, _bc, _eos, 
       force_, time_, dt);
    _Vertices = res.edges;
    _Cells = res.cells;
    _ConservedIntensive = res.intensive;
    _ConservedExtensive = res.extensive;
  }
  else
    throw UniversalError("Unsupported Runge Kutta order");

  time_ += dt;
  cycle_++;
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
  class EntropyCalculator: public Index2Member<double>
  {
  public:

    EntropyCalculator(const vector<double>& grid,
		      const vector<Primitive>& cells,
		      const EquationOfState& eos):
      grid_(grid), cells_(cells), eos_(eos) {}

    size_t getLength(void) const
    {
      return cells_.size();
    }

    double operator()(size_t i) const
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

  class PrimitiveRetriever: public Index2Member<Primitive>
  {
  public:

    PrimitiveRetriever(const vector<Conserved>& intensive,
		       const vector<Conserved>& extensive,
		       const EquationOfState& eos,
		       const vector<double> entropies,
		       double threshold,
		       bool active):
      intensive_(intensive),
      extensive_(extensive),
      eos_(eos),
      entropies_(entropies),
      threshold_(threshold),
      active_(active) {}

    size_t getLength(void) const
    {
      return intensive_.size();
    }

    Primitive operator()(size_t i) const
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


void hdsim1D::enableColdFlows(double thres)
{
  cold_flows_.activate(thres);
}
