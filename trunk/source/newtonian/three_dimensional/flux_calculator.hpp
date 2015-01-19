/*! \file flux_calculator.hpp
  \brief Abstract class for flux calculator
  \author Almog Yalinewich
 */

#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include "conserved_3d.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"

//! \brief Abstract class for flux calculator
class FluxCalculator
{
public:

  /*! \brief Calculates the fluxes
    \param tess Tessellation
    \param cells Computational cells
    \param eos Equation of state
    \param point_velocities Velocities of the mesh generating point
    \return Fluxes
   */
  virtual vector<Conserved3D> operator()
  (const Tessellation3D& tess,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const vector<Vector3D>& point_velocities) const = 0;

  //! \brief Class destructor
  virtual ~FluxCalculator(void);
};

#endif // FLUX_CALCULATOR_HPP
