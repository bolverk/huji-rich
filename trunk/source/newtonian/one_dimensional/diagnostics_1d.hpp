/*! \file diagnostics_1d.hpp
  \brief Output method for the 1d simulation
  \author Almog Yalinewich
  \deprecated We now use hdf5 as output format, so this file should be removed
 */

#ifndef DIAGNOSTICS_1D
#define DIAGNOSTICS_1D 1

#include <iostream>
#include <string>
#include <fstream>
#include "hdsim.hpp"

using namespace std;

/*! \brief Retrieves a property of the computational cell
  \param sim Hydrodynamic simulation
  \param i Cell index
  \param property Name of the property
  \return Cell property
 */
double cell_property(hdsim1D const& sim,
		     int i,
		     string const& property);

/*! \brief Creates a list of the same property from all cells
  \param sim Hydrodynamic simulation
  \param property Name of the property
  \return List of numbers
 */
vector<double> cells_property(hdsim1D const& sim,
			      string const& property);

/*! \brief Writes cells' properties to a file
  \param sim Hydrodynamic simulation
  \param property Name of property
  \param fname Name of output file
  \param pres Precision (number of digits)
 */
void write_cells_property(hdsim1D const& sim,
			  string const& property,
			  string const& fname,
			  int pres = 6);

//! \brief Class for managing actions to work at regular intervals of (virtual) time
class Timer
{
public:

  /*! \brief Class constructor
    \param t_next Next time for timer to go off
    \param dt Difference between consecutive time
   */
  Timer(double t_next, double dt);

  /*! \brief Check whether it is time
    \param t Current time
    \return True if it is time, false otherwise
   */
  bool isTime(double t);

  //! \brief Return the number of times the Timer went off
  int getCycle(void) const;

private:
  double t_next_;
  const double dt_;
  int cycle_;
};

#endif // DIAGNOSTICS_1D
