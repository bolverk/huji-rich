/*! \file cell_updater_1d.hpp
  \author Almog Yalinewich
  \brief Base class for cell update schemes
*/

#ifndef CELL_UPDATER_1D_HPP
#define CELL_UPDATER_1D_HPP 1

#include <vector>
#include "source/newtonian/two_dimensional/computational_cell_2d.hpp"
#include "source/newtonian/common/hydrodynamic_variables.hpp"
#include "source/newtonian/common/equation_of_state.hpp"
#include "physical_geometry_1d.hpp"
#include "../two_dimensional/extensive.hpp"
#include "simulation_state_1d.hpp"

using std::vector;

//! \brief Abstract class for the cell update scheme
class CellUpdater1D
{
public:
  /*! \brief Updates cells
    \param pg Physical geometry
    \param extensives Extensive variables
    \param old Old cells
    \param eos Equation of state
    \return New computational cells
  */
  virtual vector<ComputationalCell> operator()
  (const PhysicalGeometry1D& pg,
   const vector<Extensive>& extensives,
   const SimulationState1D& old,
   const EquationOfState& eos) const = 0;
    
  virtual ~CellUpdater1D(void);
};

//! \brief A simple cell updater
class SimpleCellUpdater1D: public CellUpdater1D
{
public:

  //! \brief Class constructor
  SimpleCellUpdater1D(void);
        
  vector<ComputationalCell> operator()
  (const PhysicalGeometry1D& pg,
   const vector<Extensive>& extensives,
   const SimulationState1D& old,
   const EquationOfState& eos) const;
         
  ~SimpleCellUpdater1D(void);
         
};

#endif // CELL_UPDATER_1D_HPP
