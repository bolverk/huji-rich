/*! \file main_loop_1d.hpp
  \brief Standard simulation time advance loop
  \author Almog Yalinewich
 */

#ifndef MAIN_LOOP_HPP
#define MAIN_LOOP_HPP 1

#include <string>
#include "../one_dimensional/hdsim.hpp"

using std::string;

//! \brief Functions and classes for managing simulation time advance loops
namespace simulation1d{

  //! \brief Abstract type for termination condition
  class TerminationCondition
  {
  public:

    /*! \brief Determines if a simulation should continue running
      \param sim Hydrodynamic simulation
      \return True if the simulation should run more
     */
    virtual bool operator()(hdsim1D const& sim) = 0;

    virtual ~TerminationCondition(void);

  private:
  };

  //! \brief Terminates the simulation after a certain time is reached
  class SafeTimeTermination: public TerminationCondition
  {
  public:

    /*! \brief Class constructor
      \param termination_time Time to stop simulation
      \param max_cycles Upper bound on the number of time advance cycles
     */
    SafeTimeTermination(double termination_time,
			int max_cycles);

    bool operator()(hdsim1D const&  sim);

  private:
    
    const double termination_time_;
    const int max_cycles_;
  };

  //! \brief Abstract class for diagnostics function
  class DiagnosticsFunction
  {
  public:

    /*! \brief Perform diagnostics
      \param sim Hydrodynamic simulation
     */
    virtual void diagnose(hdsim1D const& sim) = 0;
    
    virtual ~DiagnosticsFunction(void);
  };

  //! \brief Writes the time to a file after each time advance cycle
  class WriteTime: public DiagnosticsFunction
  {
  public:

    /*! \brief Class constructor
      \param fname Output file name
     */
    explicit WriteTime(string const& fname);

    /*! \brief Performs diagnostics
      \param sim Hydrodynamic simulation
     */
    void diagnose(hdsim1D const& sim);

  private:
    const string fname_;
  };

  /*! \brief Main simulation time advance loop
    \param sim Hydrodynamic simulation
    \param term_cond Termination condition
    \param time_order Time integration order
    \param diag Diagnostic function
   */
  void main_loop(hdsim1D& sim,
		 TerminationCondition& term_cond,
		 int time_order,
		 DiagnosticsFunction* diag = 0);

  /*! \brief Main simulation time advance loop
    \param sim Hydrodynamic simulation
    \param final_time Final simulation time
    \param max_iter Maximum number of time advance cycles
    \param time_order Time integration order
    \param time_log Name of output file
   */
  void main_loop(hdsim1D& sim,
		 double final_time,
		 int max_iter=1e6,
		 int time_order=1,
		 string const& time_log="");
}

#endif // MAIN_LOOP_HPP
