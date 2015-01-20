/*! \file hdf5_logger.hpp
  \brief Method for dumping tessellation data to hdf5 file
 */

#ifndef HDF5_LOGGER_HPP
#define HDF5_LOGGER_HPP 1

#include "voronoi_logger.hpp"

using std::string;

class HDF5Logger: public voronoi_loggers::VoronoiLogger
{
public:

  HDF5Logger(const string& fname);

  void output(const VoronoiMesh& v);

  void output(const Tessellation& v);

private:
  const string fname_;
};

#endif // HDF5_LOGGER_HPP

