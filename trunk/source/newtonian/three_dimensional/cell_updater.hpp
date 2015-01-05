#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"

class CellUpdater
{
public:

  virtual ComputationalCell operator()
  (const Conserved3D& intensive,
   const EquationOfState& eos) const = 0;

  virtual ~CellUpdater(void);
};

#endif // CELL_UPDATER_HPP
