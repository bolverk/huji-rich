#include "computational_cell.hpp"

ComputationalCell::ComputationalCell(void):
  density(0),
  pressure(0),
  velocity(0,0),
  tracers(),
  stickers() {}
