#include "condition_action_sequence.hpp"

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
