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
#include "../../misc/utils.hpp"
#include <cmath>

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
void BinOutput(string location,
	       hdsim const& sim,
	       Tessellation const& V,
	       bool floatprecision=true);

/*! \brief Calculates the total extensive conserved variables of the entire computational domain
  \param sim Hydrodynamic simulation
  \return The total conserved of the simulation
*/
Conserved total_conserved(hdsim const& sim);

/*! \brief Displays the UniversalError information
  \param eo Th UniversalError
*/
void DisplayError(UniversalError const& eo);

/*! \brief Returns the ConvexHull for a set of points
  \param result The set of convex hull points
  \param tess The tessellation
  \param index The index of the cell for which to calculate the convex hull
*/
void ConvexHull(vector<Vector2D> &result,Tessellation const* tess,int index);

#endif // DIAGNOSTICS_HPP
