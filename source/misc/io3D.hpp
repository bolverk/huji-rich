/*! \file io3D.hpp
\brief A collection of simple input / output methods for Vector3D
\author Elad Steinberg
*/

#ifndef IO3D_HPP
#define IO3D_HPP 1

#include "simple_io.hpp"
#include "../3D/GeometryCommon/Vector3D.hpp"

/*! \brief Writes a vector of Vector3D to a file
\param vec The vector to write
\param fname Name of the file
*/
void write_vec3d(std::vector<Vector3D> const&vec, std::string const& fname);

/*! \brief Reads a  a vector of Vector3D from a file
\param fname Name of the file
\return The vector
*/
std::vector<Vector3D> read_vec3d(std::string fname);
#endif // IO3D_HPP
