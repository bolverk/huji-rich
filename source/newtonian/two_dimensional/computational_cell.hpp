#ifndef COMPUTATIONAL_CELL_HPP
#define COMPUTATIONAL_CELL_HPP 1

#include <map>
#include <string>
#include "../../tessellation/geometry.hpp"

class ComputationalCell
{
public:

  double density;
  double pressure;
  Vector2D velocity;

  std::map<std::string,double> tracers;

  std::map<std::string,bool> stickers;
};

#endif // COMPUTATIONAL_CELL_HPP
