/*! \file BinaryLogger.hpp
\brief Class for a logger for the tessellation in 3D that outputs a binary file. The output is as follows:
1) size_t Number of points (N)
2) For each point x,y,z coordinates (3 doubles)
3) size_t Number of faces (Nfaces)
4) For each face, the neighbors (2*size_t), the number of vertices (size_t) and the vertices (3*double*number of vertices)
5) The list of faces per cell, the number of neighbors (size_t) and then the indeces of the faces (size_t*number of neighbors)
\author Elad Steinberg
*/

#ifndef BINARYLOGGER_HPP
#define BINARYLOGGER_HPP 1

#include "TessLogger.hpp"
#include <fstream>
#include <string>
#include "../../misc/universal_error.hpp"
#include "../../misc/simple_io.hpp"

//! \brief Tessellation debugging diagnostic that writes the data to a binary file
class BinaryLogger
{
private:
	string filename_;
public:
	/*!
	\brief Class constructor
	\param filename The name of the output file
	*/
	BinaryLogger(string const& filename);

	//! \brief class destructor
	~BinaryLogger(void);

	/*!
	\brief Creates a log file of the tessellation
	\param tess The tessellation to log
	*/
	void Log(Tessellation3D const& tess)const;
};
#endif //BINARYLOGGER_HPP
