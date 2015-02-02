/*! \file consecutive_snapshots.hpp
  \brief A diagnostics class that writes snapshots at regular intervals
  \author Almog Yalinewich
 */

#ifndef CONSECUTIVE_SNAPSHOTS_HPP
#define CONSECUTIVE_SNAPSHOTS_HPP 1

#include <memory>
#include "main_loop_2d.hpp"
#include "index2filename.hpp"

//! \brief A diagnostic class that writes snapshots at regular intervals
class ConsecutiveSnapshots: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param dt Time step
    \param init_time Initial time
    \param counter Counts how many times time advance was called
   */
  ConsecutiveSnapshots(double dt,
		       Index2FileName* i2fn,
		       double init_time=0,
		       int counter=0);

  void operator()(hdsim const& sim);

private:
  double next_time_;
  const double dt_;
  int counter_;
  std::auto_ptr<Index2FileName> i2fn_;
};

#endif // CONSECUTIVE_SNAPSHOTS_HPP
