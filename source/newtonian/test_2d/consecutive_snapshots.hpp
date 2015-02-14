/*! \file consecutive_snapshots.hpp
  \brief A diagnostics class that writes snapshots at regular intervals
  \author Almog Yalinewich
 */

#ifndef CONSECUTIVE_SNAPSHOTS_HPP
#define CONSECUTIVE_SNAPSHOTS_HPP 1

#include <memory>
#include "main_loop_2d.hpp"
#include "index2filename.hpp"
#include "trigger.hpp"

//! \brief A diagnostic class that writes snapshots at regular intervals
class ConsecutiveSnapshots: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param dt Time step
    \param init_time Initial time
    \param counter Counts how many times time advance was called
   */
  ConsecutiveSnapshots(Trigger* trigger,
		       Index2FileName* i2fn);

  void operator()(hdsim const& sim);

private:
  std::auto_ptr<Trigger> trigger_;
  std::auto_ptr<Index2FileName> i2fn_;
  int counter_;
};

#endif // CONSECUTIVE_SNAPSHOTS_HPP
