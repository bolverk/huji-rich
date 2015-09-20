/*! \file simple_flux_calculator.hpp
  \author Almog Yalinewich
  \brief Simple flux calculator
 */

#ifndef PERIODIC_BC_HPP
#define PERIODIC_BC_HPP 1

#include "flux_calculator_2d.hpp"
#include "../common/riemann_solver.hpp"

//! \brief Simple flux calculator
class PeriodicBC: public FluxCalculator
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
   */
  explicit PeriodicBC
  (const RiemannSolver& rs);

  vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double time,
   const double dt) const;

private:
  const RiemannSolver& rs_;

  Conserved calcHydroFlux
  (const Tessellation& tess,
   const vector<Vector2D>& point_velocities,
   const vector<ComputationalCell>& cells,
   const EquationOfState& eos,
   const size_t i) const;
};

#endif // SIMPLE_FLUX_CALCULATOR_HPP
