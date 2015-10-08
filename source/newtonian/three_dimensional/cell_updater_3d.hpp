/*! \file cell_updater.hpp
  \brief Abstract class for cell update scheme
  \author Almog Yalinewich
 */

#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"

using three_dimenssional::ComputationalCell;

//! \brief Abstract clas for cell update scheme
class CellUpdater
{
public:

  /*! \brief Calculates the computational cell
    \param intensive Intensive conserved variables (per volume)
    \param eos Equation of state
    \return Computational cell
   */
  virtual ComputationalCell operator()
  (const Conserved3D& intensive,
   const EquationOfState& eos) const = 0;

  //! \brief Class destructor
  virtual ~CellUpdater(void);
};

#endif // CELL_UPDATER_HPP
