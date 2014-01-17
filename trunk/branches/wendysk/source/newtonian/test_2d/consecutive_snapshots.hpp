#ifndef CONSECUTIVE_SNAPSHOTS_HPP
#define CONSECUTIVE_SNAPSHOTS_HPP 1

#include "main_loop_2d.hpp"

class ConsecutiveSnapshots: public DiagnosticFunction
{
public:

  ConsecutiveSnapshots(double dt);

  void diagnose(hdsim const& sim);

private:
  double next_time_;
  const double dt_;
  int counter_;
};

#endif // CONSECUTIVE_SNAPSHOTS_HPP
