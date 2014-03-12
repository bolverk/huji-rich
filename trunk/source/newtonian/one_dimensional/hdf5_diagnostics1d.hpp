/*! \file hdf5_diagnostics1d.hpp
  \brief Output method formatted in hdf5
  \author Almog Yalinewich
 */ 

#include <H5Cpp.h>
#include <string>
#include "hdsim.hpp"

using std::string;

//! \brief Diagnostics for one dimensional simulation
namespace diagnostics1d{

  /*! \brief Writes all hydrodynamic data to an hdf5 file
    \param sim Hydrodynamic simulation
    \param fname File name
   */
  void write_snapshot_to_hdf5(hdsim1D const& sim,
			      string const& fname);
}
