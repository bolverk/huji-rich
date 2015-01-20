/*! \file hdf5_utils.hpp
  \brief Higher level hdf5 utilities
  \author Almog Yalinewich
 */

#ifndef HDF5_UTILS_HPP
#define HDF5_UTILS_HPP 1

#include <string>
#include <vector>

using std::string;
using std::vector;
using std::pair;

//! \brief Facilitates writing hdf5 files
class HDF5Shortcut
{
public:

  /*! \brief Class constructor
    \param fname Name of hdf5 file
   */
  HDF5Shortcut(const string& fname);

  /*! \brief adds dataset
    \param field_name Name of dataset
    \param data_set Array of data
    \return Self reference
   */
  HDF5Shortcut& operator()(const string& field_name,
			   const vector<double>& data_set);
  
  /*! \brief adds dataset
    \param field_name Name of dataset
    \param data_set Array of data
    \return Self reference
   */
  HDF5Shortcut& operator()(const string& field_name,
			   const vector<int>& data_set);

  //! \brief Class destructor. This is the stage when the file is written
  ~HDF5Shortcut(void);

private:
  const string fname_;
  vector<pair<string,vector<double> > > double_data_;
  vector<pair<string,vector<int> > > int_data_;
};

#endif // HDF5_UTILS_HPP
