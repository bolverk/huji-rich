/*! \file flux_calculator_3d.hpp
  \brief Abstract class for flux calculator
  \author Almog Yalinewich
 */

#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include "conserved_3d.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"
#include "../common/equation_of_state.hpp"
#include "RiemannSolver3D.hpp"

//! \brief Abstract class for flux calculator
class FluxCalculator3D
{
public:

  /*! \brief Calculates the fluxes
  \param extensives The extensive variables
    \param tess Tessellation
    \param cells Computational cells
    \param eos Equation of state
    \param edge_velocities Velocities of the edges
    \param fluxes THe fluxes given as output
	\param dt The timestep
	\param time The time
	\param tracerstickernames The names of the stickers and tracers
   */
	virtual std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > operator()(vector<Conserved3D>& fluxes, const Tessellation3D& tess, const vector<Vector3D>& edge_velocities,
	  const vector<ComputationalCell3D>& cells,const vector<Conserved3D>& extensives,const EquationOfState& eos,
	  const double time, const double dt,TracerStickerNames const& tracerstickernames) const = 0;

  //! \brief Class destructor
  virtual ~FluxCalculator3D(void);
};

void RotateSolveBack3D(Vector3D const& normal, ComputationalCell3D const& left, ComputationalCell3D const& right,
	Vector3D const& face_velocity,RiemannSolver3D const& rs, Conserved3D &res,EquationOfState const& eos,TracerStickerNames const& tsn);

#endif // FLUX_CALCULATOR_HPP
