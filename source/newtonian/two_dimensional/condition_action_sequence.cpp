#include "condition_action_sequence.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"

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
     const Vector2D& edge_velocity,
     const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence)
  {
    for(size_t i=0;i<sequence.size();++i){
      const pair<bool,bool> flag_aux = (*sequence[i].first)
	(edge,tess,cells);
      if(flag_aux.first)
	return (*sequence[i].second)
	  (edge,tess, edge_velocity,cells,eos,flag_aux.second);
    }
    throw UniversalError("Error in ConditionActionSequence");
  }
}

vector<Extensive> ConditionActionSequence::operator()
(const Tessellation& tess,
 const vector<Vector2D>& edge_velocities,
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
		  edge_velocities[i],
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
}

Extensive RegularFlux::operator()
(const Edge& edge,
 const Tessellation& tess,
	const Vector2D& edge_velocity,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const bool /*aux*/) const
{
	assert(edge.neighbors.first >= 0 && tess.GetOriginalIndex(edge.neighbors.first) !=
		tess.GetOriginalIndex(edge.neighbors.second) && edge.neighbors.second >= 0);
  const Vector2D p = normalize
    (edge.vertices.second -
     edge.vertices.first);
  const Vector2D n = normalize
    (tess.GetMeshPoint(edge.neighbors.second)-
     tess.GetMeshPoint(edge.neighbors.first));
  const double v =
    ScalarProd(n,edge_velocity);
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

RigidWallFlux::RigidWallFlux
(const RiemannSolver& rs):
  rs_(rs) {}

namespace{
  pair<Primitive,Primitive> rigid_wall_states
    (const Primitive& state,
     const Vector2D& p,
     const bool aux)
  {
    if(aux){
      const Primitive left = state;
      const Primitive right = reflect(left,p);
      return pair<Primitive,Primitive>(left,right);
    }
    else{
      const Primitive right = state;
      const Primitive left = reflect(right,p);
      return pair<Primitive,Primitive>(left,right);
    }
  }
}

Extensive RigidWallFlux::operator()
(const Edge& edge,
 const Tessellation& tess,
	const Vector2D& /*edge_velocity*/,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const bool aux) const
{
  if(aux)
    assert(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo());
  else
    assert(edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
  const Vector2D p = normalize
    (edge.vertices.second - edge.vertices.first);
  const Vector2D n = 
    normalize
    (remove_parallel_component
     (aux ?
      edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first) :
      tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
      p));
  const double v = 0;
  const pair<Primitive,Primitive> left_right =
    rigid_wall_states
    (convert_to_primitive
     (cells.at
      (static_cast<size_t>
       (aux ? edge.neighbors.first  : edge.neighbors.second)),
      eos),
     p, aux);
  const Conserved c = rotate_solve_rotate_back
    (rs_,
     left_right.first,
     left_right.second,
     v,n,p);
  return conserved_to_extensive
    (c,
     cells.at(static_cast<size_t>(aux ? edge.neighbors.first : edge.neighbors.second)));
}

FreeFlowFlux::FreeFlowFlux(const RiemannSolver& rs):
  rs_(rs) {}

Extensive FreeFlowFlux::operator()
(const Edge& edge,
 const Tessellation& tess,
	const Vector2D& /*edge_velocity*/,
 const vector<ComputationalCell>& cells,
 const EquationOfState& eos,
 const bool aux) const
{
#ifndef RICH_MPI
  if(aux)
    assert(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo());
  else
    assert(edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
#endif //RICH_MPI
  const Vector2D p = normalize
    (edge.vertices.second - edge.vertices.first);
  const Vector2D n = 
    normalize
    (remove_parallel_component
     (aux ?
      edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first) :
      tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
      p));
  const double v = 0;
  const Primitive state =
    convert_to_primitive
    (cells.at
     (static_cast<size_t>
      (aux ? edge.neighbors.first : edge.neighbors.second)),
     eos);
  const Conserved c = rotate_solve_rotate_back
    (rs_,
     state, state,
     v,n,p);
  return conserved_to_extensive
    (c,
     cells.at(static_cast<size_t>(aux ? edge.neighbors.first : edge.neighbors.second)));
}

IsBoundaryEdge::IsBoundaryEdge(void) {}

pair<bool,bool> IsBoundaryEdge::operator()
(const Edge& edge,
 const Tessellation& tess,
 const vector<ComputationalCell>& /*cells*/) const
{
#ifdef RICH_MPI
	if(tess.GetOriginalIndex(edge.neighbors.first)!=tess.GetOriginalIndex(edge.neighbors.second))
		return pair<bool, bool>(false, false);
	else
	{
		if (edge.neighbors.first < tess.GetPointNo())
			return pair<bool, bool>(true, true);
		else
			return pair<bool, bool>(true, false);
	}
#endif
  if(edge.neighbors.first<0 || edge.neighbors.first>=tess.GetPointNo()){
    assert(edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo());
    return pair<bool,bool>(true,false);
  }
  if(edge.neighbors.second<0 || edge.neighbors.second>=tess.GetPointNo())
    return pair<bool,bool>(true,true);
  return pair<bool,bool>(false,false);
}

IsBulkEdge::IsBulkEdge(void) {}

pair<bool,bool> IsBulkEdge::operator()
(const Edge& edge,
 const Tessellation& tess,
 const vector<ComputationalCell>& /*cells*/) const
{
  return pair<bool,bool>
    (edge.neighbors.first >= 0 &&
     edge.neighbors.second >= 0 &&
     ((edge.neighbors.first < tess.GetPointNo() &&
		 edge.neighbors.second < tess.GetPointNo() ) || 
		(tess.GetOriginalIndex(edge.neighbors.first)!=
		 tess.GetOriginalIndex(edge.neighbors.second))),
     false);
}

RegularSpecialEdge::RegularSpecialEdge(const string& sticker_name):
  sticker_name_(sticker_name) {}

pair<bool,bool> RegularSpecialEdge::operator()
(const Edge& edge,
 const Tessellation& /*tess*/,
 const vector<ComputationalCell>& cells) const
{
   if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.first)).stickers,sticker_name_)){
    if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,sticker_name_))
      return pair<bool,bool>(false,false);
    return pair<bool,bool>(true,false);
  }
  if(safe_retrieve(cells.at(static_cast<size_t>(edge.neighbors.second)).stickers,sticker_name_))
    return pair<bool,bool>(true,true);
  return pair<bool,bool>(false,false);
}
