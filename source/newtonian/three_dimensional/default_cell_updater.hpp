/*! \file default_cell_updater.hpp
  \brief Default cell update scheme
  \author Almog Yalinewich
 */

#ifndef DEFAULT_CELL_UPDATER_HPP
#define DEFAULT_CELL_UPDATER_HPP 1

#include "cell_updater.hpp"

//! \brief Default scheme for cell update
class DefaultCellUpdater: public CellUpdater
{
public:

  //! \brief Class constructor
  DefaultCellUpdater(void);

  ComputationalCell operator()
  (const Conserved3D& intensive,
   const EquationOfState& eos) const;
};

#endif // DEFAULT_CELL_UPDATER_HPP
