/*! \file hdf5_diagnostics.hpp
  \brief Simulation output to hdf5 file format
  \author Elad Steinberg
 */

#ifndef HDF5_DIAG
#define HDF5_DIAG 1

#include <H5Cpp.h>
#include <string>
#include "hdsim2d.hpp"
#include "../../misc/int2str.hpp"
#include "diagnostics.hpp"
#include "Reset.hpp"

//! \brief High level function for writing hdf5 files
class HDF5Shortcut
{
public:

  /*! \brief Class constructor
    \param fname Name of hdf5 file
   */
  HDF5Shortcut(const string& fname);

  /*! \brief Adds data set
    \param field_name Name of dataset
    \param array Data
   */
  void operator()(const string& field_name,
		  const vector<double>& array);

  //! \brief Class destructor. Dumps data to file
  ~HDF5Shortcut(void);

private:
  const string fname_;
  vector<std::pair<string, vector<double> > > data_;
};

/*!
\brief Writes the simulation data into an HDF5 file
\param sim The hdsim class of the simulation
\param fname The name of the output file
*/
void write_snapshot_to_hdf5(hdsim const& sim,string const& fname);
/*!
\brief Reads an HDF5 snapshot file in order to restart the simulation
\param dump The dump data structure, should be when passed
\param fname The path to the HDF5 file
\param eos The equation of state
*/
void read_hdf5_snapshot(ResetDump &dump,string const& fname,EquationOfState
	const* eos);
/*!
\brief Converts an HDF5 snapshot file to the RICH custom reset binary format
\param input The path to the HDF5 file
\param output The path to the new binary file
*/
void ConvertHDF5toBinary(string const& input, string const& output);

#endif // HDF5_DIAG
