/*! \file Reset.hpp
  \brief Functions and classes for restarting a simulation from an output files
  \author Elad Steinberg
 */

#ifndef RESET_HPP
#define RESET_HPP 1
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"
#include <fstream>
#include <string>

/*!
\brief Outputs the simulation data into a restart file
\param location The filename of the dump
\param sim The sim data
*/
void ResetOutput(string location,hdsim const& sim);
/*!
\brief Reads the simulation data from a restart file
\param location The filename of the dump
\param dump The dump class that is overwritten
\param eos The equation of state
*/
void ResetRead(string location,ResetDump &dump,EquationOfState const* eos);

#endif // RESET_HPP
