/*! \file BinaryLogger.hpp
\brief Class for a logger for the tessellation in 3D that outputs a binary file
\author Elad Steinberg
*/

#ifndef BINARYLOGGER_HPP
#define BINARYLOGGER_HPP 1

#include "TessLogger.hpp"
#include <fstream>
#include "../../misc/universal_error.hpp"
#include "../../misc/simple_io.hpp"

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

	void Log(Tessellation3D const& tess)const;
}
#endif //BINARYLOGGER_HPP
