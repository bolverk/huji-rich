#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include "computational_cell.hpp"

class CellUpdater
{
public:

  virtual ComputationalCell operator()
  (const Conserved3D& intensive,
   const EquationOfState& eos) const = 0;
};

#endif // CELL_UPDATER_HPP
