/*! \file computational_cell.hpp
  \author Almog Yalinewich
  \brief Computational cell
 */

#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

#include <map>
#include <string>
#include "../../tessellation/geometry.hpp"

//! \brief Computational cell
class ComputationalCell
{
public:

  //! \brief Density
  double density;

  //! \brief Pressure
  double pressure;

  //! \brief Velocity
  Vector2D velocity;

  //! \brief Tracers (can transfer from one cell to another)
  std::map<std::string,double> tracers;

  //! \brief Stickers (stick to the same cell)
  std::map<std::string,bool> stickers;
};

#endif // COMPUTATIONAL_CELL_HPP
