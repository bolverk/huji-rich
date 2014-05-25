/*! \file diagnostics.hpp
  \author Almog Yalinewich
  \brief Contains function that write simulation data to a file
 */

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP 1

#include <iostream>
#include <fstream>
#include <string>
#include "hdsim2d.hpp"

using namespace std;

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
void BinOutput(string location,hdsim const& sim,Tessellation const& V,bool
	floatprecision=true);

/*! \brief Goes through all cells and extracts one of their properties to a list
  \param sim Hydrodynamic simulation
  \param property Name of the property
  \return The extracted property
 */
vector<double> cells_property(hdsim const& sim,
			      string const& property);

/*! \brief Cycles through all cells and writes out one of their properties
  \param sim Hydrodynamic simulation
  \param property Name of the property to be written
  \param fname Name of output file
 */
void write_cells_property(hdsim const& sim,
			  string const& property,
			  string const& fname);

/*! \brief Replaces all occurances of a character inside a string with a new characeter
  \param base Original string
  \param char_old Character to be replaced
  \param char_new New character
  \return The new string
 */
string replace_all(string const& base,
		   char char_old, char char_new);

/*! \brief Goes through all the cell sides, writes out their vertices and neighbors
  \param sim Hydrodynamic simulation
  \param fname Name of output file
 */
void write_edges_and_neighbors(hdsim const& sim,
			       string const& fname);

/*! \brief Writes the position of each mesh generating point
  \param sim Hydrodynamic simulation
  \param fname Name of the output file
 */
void write_generating_points(hdsim const& sim,
			     string const& fname);

/*! \brief Writes a two dimensional array of numbers to a file
  \param data Two dimensional array
  \param fname Output file name
 */
void write_array_2d(vector<vector<double> > const& data,
		    string const& fname);

/*! \brief Calculates the total extensive conserved variables of the entire computational domain
  \param sim Hydrodynamic simulation
  \return The total conserved of the simulation
 */
Conserved total_conserved(hdsim const& sim);

/*! \brief Displays the UniversalError information
	\param eo Th UniversalError
	\param cycle_number The simulation iteration number
*/

void DisplayError(UniversalError const& eo,int cycle_number);

/*! \brief Returns the ConvexHull for a set of points
\param result The set of convex hull points
\param tess The tessellation
\param index The index of the cell for which to calculate the convex hull
*/
void ConvexHull(vector<Vector2D> &result,Tessellation const* tess,int index);

#endif // DIAGNOSTICS_HPP
