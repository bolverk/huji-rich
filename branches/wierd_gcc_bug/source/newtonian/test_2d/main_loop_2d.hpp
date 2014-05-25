#ifndef MAIN_LOOP_2D_HPP
#define MAIN_LOOP_2D_HPP 1

#include <string>
#include "../two_dimensional/hdsim2d.hpp"

using namespace std;

class DiagnosticFunction
{
public:

  virtual void diagnose(hdsim const& sim) = 0;

  virtual ~DiagnosticFunction(void);
};

class WriteTime: public DiagnosticFunction
{
public:

  WriteTime(string const& fname);

  void diagnose(hdsim const& sim);

private:
  string fname_;
};

class TerminationCondition
{
public:

  virtual bool should_continue(hdsim const& sim) = 0;

  virtual ~TerminationCondition(void) = 0;
};

class SafeTimeTermination: public TerminationCondition
{
public:

  SafeTimeTermination(double termination_time,
		      int max_cycles);

  bool should_continue(hdsim const& sim);

private:
  const double termination_time_;
  const int max_cycles_;
};

class CycleTermination: public TerminationCondition
{
public:

  CycleTermination(int max_cycles);

  bool should_continue(hdsim const& sim);

private:

  const int max_cycles_;
};

namespace simulation2d{
  void main_loop(hdsim& sim,
		 TerminationCondition& term_cond,
		 int time_order = 1,
		 DiagnosticFunction* diagfunc = 0);
}

#endif // MAIN_LOOP_2D_HPP
