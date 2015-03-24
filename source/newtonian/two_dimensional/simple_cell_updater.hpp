/*! \file simple_cell_updater.hpp
  \author Almog Yalinewich
  \brief Simple cell updater
 */

#ifndef SIMPLE_CELL_UPDATER_HPP
#define SIMPLE_CELL_UPDATER_HPP 1

#include "cell_updater.hpp"
#include "../../misc/utils.hpp"

using std::vector;

//! \brief Simple cell updater
class SimpleCellUpdater: public CellUpdater
{
public:

SimpleCellUpdater(void);

vector<ComputationalCell> operator()
(const Tessellation& tess,
   const PhysicalGeometry& pg,
   const EquationOfState& eos,
   const vector<Extensive>& extensives,
   const vector<ComputationalCell>& old) const;
};

#endif // SIMPLE_CELL_UPDATER_HPP
