#ifndef CONDITION_ACTION_SEQUENCE_HPP
#define CONDITION_ACTION_SEQUENCE_HPP 1

#include "flux_calculator_2d.hpp"
#include "../common/riemann_solver.hpp"

class ConditionActionSequence: public FluxCalculator
{
public:

  class Condition
  {
  public:

    virtual pair<bool, bool> operator()
    (const Edge& edge,
     const Tessellation& tess,
     const vector<ComputationalCell>& cells) const = 0;

    virtual ~Condition(void);
  };

  class Action
  {
  public:

    virtual Extensive operator()
    (const Edge& edge,
     const Tessellation& tess,
     const vector<Vector2D>& point_velocities,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const bool aux) const = 0;

    virtual ~Action(void);
  };

  ConditionActionSequence
  (const vector<pair<const Condition*, const Action*> >& sequence);

  ~ConditionActionSequence(void);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double time,
   const double dt) const;

private:
   const vector<pair<const Condition*, const Action*> > sequence_;
};

class RegularFlux: public ConditionActionSequence::Action
{
public:

  RegularFlux(const RiemannSolver& rs);

  Extensive operator()
  (const Edge& edge,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const bool aux) const;

private:

  const RiemannSolver& rs_;
};

class RigidWallFlux: public ConditionActionSequence::Action
{
public:

  RigidWallFlux(const RiemannSolver& rs);

  Extensive operator()
  (const Edge& edge,
   const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const bool aux) const;

private:
  const RiemannSolver& rs_;
};

#endif // CONDITION_ACTION_SEQUENCE_HPP
