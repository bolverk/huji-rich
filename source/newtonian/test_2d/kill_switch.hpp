/*! \file kill_switch.hpp
  \author Almog Yalinewich
  \brief Allows user to manually terminate the main loop
 */

#ifndef KILL_SWITCH_HPP
#define KILL_SWITCH_HPP 1

#include "main_loop_2d.hpp"

//! \brief Wrapper around another termination condition that allows the user to terminate the main loop manually
class KillSwitch: public TerminationCondition
{
public:

  /*! \brief Class constructor
    \param fname Kill switch file. When the program first starts, a file with this name will be created and the value 0 be written to it. If, at some point, the value in this file will be different from 0, then the main loop should stop.
    \param tc Original termination condition
   */
  KillSwitch(const string& fname,
	     TerminationCondition& tc);

  bool operator()(const hdsim& sim);

private:

  const string fname_;
  TerminationCondition& tc_;
};

#endif // KILL_SWITCH_HPP
