/*! \file diagnostics.hpp
  \author Almog Yalinewich
  \brief Contains function that write simulation data to a file
*/

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP 1

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "hdsim2d.hpp"
#include "../../misc/utils.hpp"
#include "../../tessellation/ConvexHull.hpp"
#include "physical_geometry.hpp"

using std::string;

/*!
  \brief Outputs the simulation data into a binary file (float precision) as follows:
  1) Int - Number of cells=n
  2) For each cell it's x and y coordinate (2 floats)*n
  3) Pressure (n floats)
  4) Density (n floats)
  5) For each cell an int=N, stating the number of vertices and then
  the vertices(2*N floats)
  6) The time - float
  7) The number of cells with tracers - int
  8) The number of tracers - int
  9) The passive tracers written cell by cell (float)
  \param location The output file's location
  \param sim The hydro sim
  \param V The tessellation
  \param floatprecision A flag to choose output in float/double precision
*/
void BinOutput(string location,
	       hdsim const& sim,
	       Tessellation const& V,
	       bool floatprecision=true);

/*! \brief Calculates the total extensive conserved variables of the entire computational domain
  \param sim Hydrodynamic simulation
  \return The total conserved of the simulation
*/
Extensive total_conserved(hdsim const& sim);

/*! \brief Calculates the total amount of tracer in computational domain
  \param sim Hydrodynamic simulation
  \param index Index of tracer
  \return Total amount of certain tracer
 */
double total_tracer(const hdsim& sim,
		    const int index);

/*! \brief Writes all the error information to a file
  \param fname File name
  \param eo Error object
 */
void write_error(const string& fname,
		 const UniversalError& eo);

/*! \brief Writes a vector of Vector2D to a binary file
\param vec The vector to write
\param filename The path to the output file
*/
void WriteVector2DToFile(vector<Vector2D> const& vec,string filename);

/*! \brief Reads a vector of Vector2D from a binary file
\return The vector that was read
\param filename The path to the output file
*/
vector<Vector2D> ReadVector2DFromFile(string filename);

#endif // DIAGNOSTICS_HPP
