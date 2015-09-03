#include "condition_action_sequence.hpp"
#include "simple_flux_calculator.hpp"

ConditionActionSequence::ConditionActionSequence
(const vector<pair<const Condition*,const Action*> >& sequence):
  sequence_(sequence) {}

ConditionActionSequence::~ConditionActionSequence(void)
{
  for(size_t i=0;i<sequence_.size();++i){
    delete sequence_[i].first;
    delete sequence_[i].second;
  }
}

namespace{
  Extensive choose_action
    (const Edge& edge,
     const Tessellation& tess,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const vector<Vector2D>& point_velocities,
     const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence)
  {
    for(size_t i=0;i<sequence.size();++i){
      const pair<bool,bool> flag_aux = (*sequence[i].first)
	(edge,tess,cells);
      if(flag_aux.first)
	return (*sequence[i].second)
	  (edge,tess,point_velocities,cells,eos,flag_aux.second);
    }
    throw "Error in ConditionActionSequence";
  }
}

vector<Extensive> ConditionActionSequence::operator()
(const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const vector<Extensive>& /*extensives*/,
 const CacheData& /*cd*/,
 const EquationOfState& eos,
 const double /*time*/,
 const double /*dt*/) const
{
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<tess.getAllEdges().size();++i)
    res[i] = choose_action
      (tess.getAllEdges()[i],
       tess,
       cells,
       eos,
       point_velocities,
       sequence_);
  return res;
}

ConditionActionSequence::Condition::~Condition(void) {}

ConditionActionSequence::Action::~Action(void) {}

RegularFlux::RegularFlux(const RiemannSolver& rs):
  rs_(rs) {}

namespace {
  Extensive conserved_to_extensive
  (const Conserved& c, const ComputationalCell& cell)
  {
    Extensive res;
    res.mass = c.Mass;
    res.momentum = c.Momentum;
    res.energy = c.Energy;
    for(boost::container::flat_map<string,double>::const_iterator it=
	  cell.tracers.begin();
	it!=cell.tracers.end();++it)
      res.tracers[it->first] = (it->second)*c.Mass;
    return res;
  }

  Vector2D normalize(const Vector2D& v)
  {
    return v/abs(v);
  }
}

Extensive RegularFlux::operator()
(const Edge& edge,
 const Tessellation& tess,
 const vector<Vector2D>& point_velocities,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const bool /*aux*/) const
{
  assert(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo());
  assert(edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
  const Vector2D p = normalize
    (edge.vertices.second -
     edge.vertices.first);
  const Vector2D n = normalize
    (tess.GetMeshPoint(edge.neighbors.second)-
     tess.GetMeshPoint(edge.neighbors.first));
  const double v =
    ScalarProd
    (n,
     tess.CalcFaceVelocity
     (point_velocities[static_cast<size_t>(edge.neighbors.first)],
      point_velocities[static_cast<size_t>(edge.neighbors.second)],
      tess.GetCellCM(edge.neighbors.first),
      tess.GetCellCM(edge.neighbors.second),
      calc_centroid(edge)));
  const Conserved c = rotate_solve_rotate_back
    (rs_,
     convert_to_primitive
     (cells.at(static_cast<size_t>(edge.neighbors.first)),eos),
     convert_to_primitive
     (cells.at(static_cast<size_t>(edge.neighbors.second)),eos),
     v,n,p);
  return conserved_to_extensive
    (c,
     c.Mass>0 ?
     cells.at(static_cast<size_t>(edge.neighbors.first)):
     cells.at(static_cast<size_t>(edge.neighbors.second)));	      
}
