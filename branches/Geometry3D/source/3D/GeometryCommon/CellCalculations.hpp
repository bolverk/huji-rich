//\file CellCalculations.hpp
//\brief Calculations for Cells
//\author Itay Zandbank

// There's no Cell class because it's too dependant of TessellationBase (requiring face indices, etc...)
// So here we just have some handle calculations

#ifndef CELL_CALCULATIONS_HPP
#define CELL_CALCULATIONS_HPP

#include <vector>

#include "Face.hpp""
#include "Tetrahedron.hpp"

//\brief Splits a cell into tetrahedra, all touching the center of the cell
std::vector<Tetrahedron> SplitCell(const std::vector<const Face *> &cell);

//\brief Calculates the volume and center-of-mass of a cell
void CalculateCellDimensions(const std::vector<const Face *> &cell, double &volume, Vector3D &centerOfMass);

#endif