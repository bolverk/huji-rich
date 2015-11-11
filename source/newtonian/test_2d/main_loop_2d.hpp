/*! \file main_loop_2d.hpp
  \brief Standard simulation time advance loop
  \author Almog Yalinewich
 */

#ifndef MAIN_LOOP_2D_HPP
#define MAIN_LOOP_2D_HPP 1

#include <string>
#include "../two_dimensional/hdsim2d.hpp"

using std::string;

//! \brief Abstract class for a diagnostic function
class DiagnosticFunction
{
public:

  /*! \brief Perform diagnostics
    \param sim Hydrodynamic simulation
   */
  virtual void operator()(const hdsim& sim) = 0;

  virtual ~DiagnosticFunction(void);
};

//! \brief Writes the time after each time step
class WriteTime: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param fname File name
   */
  explicit WriteTime(const string& fname);

  void operator()(const hdsim& sim);

private:
  string fname_;
};

//! \brief A class for writing data to a file
class WriteData : public DiagnosticFunction
{
public:

	/*! \brief Class constructor
	\param fname File name
	*/
	explicit WriteData(const string& fname);

	void operator()(const hdsim& sim);

private:
	string fname_;
};


//! \brief Abstract class for a termination condition for the main loop
class TerminationCondition
{
public:

  /*! \brief Returns true if the simulation should continue, false otherwise
    \param sim Hydrodynamic simulation
    \return True is simulation should continue running
   */
  virtual bool operator()(const hdsim& sim) = 0;

  //! \brief Virtual destructor
  virtual ~TerminationCondition(void);
};

//! \brief Terminates the simulation after a certain time was reached
class SafeTimeTermination: public TerminationCondition
{
public:

  /*! \brief Class constructor
    \param termination_time The simulation is advanced until that time is reached
    \param max_cycles Upper limit on the number of time advance cycles. An error is thrown if this number is exceeded.
   */
  SafeTimeTermination(double termination_time,
		      int max_cycles);

  bool operator()(const hdsim& sim);

private:
  const double termination_time_;
  const int max_cycles_;
};

//! \brief Terminates the simulation after a specified number of iterations
class CycleTermination: public TerminationCondition
{
public:

  /*! \brief Class constructor
    \param max_cycles Maximum number of cycles
   */
  explicit CycleTermination(int max_cycles);

  bool operator()(const hdsim& sim);

private:

  const int max_cycles_;
};

//! \brief Class for manual tweaking with the simulation data
class Manipulate
{
public:

  /*! \brief Modifies the simulation data
    \param sim Hydrodynamic simulation
   */
  virtual void operator()(hdsim& sim) = 0;

  virtual ~Manipulate(void);
};

//! \brief Functions for managing two dimensional simulations
namespace simulation2d{
  /*! \brief Simulation time advance loop
    \param sim Hydrodynamic simulation
    \param term_cond Termination condition
    \param time_advance_method Method for time advance
    \param diagfunc Diagnostic function
    \param manipulate Method for manual modification of the simulatino data
   */
  void main_loop(hdsim& sim,
		 TerminationCondition& term_cond,
		 void (hdsim::*time_advance_method)(void),
		 DiagnosticFunction* diagfunc = 0,
		 Manipulate* manipulate = 0);
}

#endif // MAIN_LOOP_2D_HPP
