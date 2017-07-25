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
#include "../two_dimensional/hdf5_diagnostics.hpp"

//! \brief A diagnostic class that writes snapshots at regular intervals
class ConsecutiveSnapshots: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param trigger Trigger function
    \param i2fn Function for choosing file names
   */
  ConsecutiveSnapshots(Trigger* trigger,
		       Index2FileName* i2fn);

    /*! \brief Class constructor
    \param trigger Trigger function
    \param i2fn Function for choosing file names
    \param appendices Additional data to be written to snapshot
   */
  ConsecutiveSnapshots(Trigger* trigger,
		       Index2FileName* i2fn,
		       const vector<DiagnosticAppendix*>& appendices);

  void operator()(const hdsim& sim);

  ~ConsecutiveSnapshots(void);

private:
  std::unique_ptr<Trigger> trigger_;
  std::unique_ptr<Index2FileName> i2fn_;
  int counter_;
  const vector<DiagnosticAppendix*> appendices_;
};

#endif // CONSECUTIVE_SNAPSHOTS_HPP
