#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include <vector>
#include "computational_cell.hpp"
#include "../../tessellation/tessellation.hpp"
#include "physical_geometry.hpp"
#include "../common/equation_of_state.hpp"
#include "extensive.hpp"

using std::vector;

class CellUpdater
{
public:

  virtual vector<ComputationalCell> operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old) const = 0;

  virtual ~CellUpdater(void);
};

#endif // CELL_UPDATER_HPP
