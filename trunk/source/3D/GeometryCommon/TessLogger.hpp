/*! \file TessLogger.hpp
\brief Abstract class for a logger for the tessellation in 3D
\author Elad Steinberg
*/

#ifndef TESSLOGGER_HPP
#define TESSLOGGER_HPP 1

#include "Tessellation3D.hpp"

//! \brief Abstract class for debug diagnostics for the tessellation
class TessLogger
{
public:
	/*!
	\brief Creates a log file of the tessellation
	\param tess The tessellation to log
	*/
	virtual void Log(Tessellation3D const& tess)const=0;

  //! \brief class destructor
  virtual ~TessLogger(void);
};
#endif //TESSLOGGER_HPP
