#ifndef MODULAR_FLUX_CALCULATOR_HPP
#define MODULAR_FLUX_CALCULATOR_HPP 1

#include "flux_calculator_2d.hpp"
#include "spatial_reconstruction.hpp"
#include "../common/riemann_solver.hpp"
#include "HydroBoundaryConditions.hpp"

//! \brief Modular flux calculator
class ModularFluxCalculator: public FluxCalculator
{
public:

  /*! \brief Class constructor
    \param sr Interpolation
    \param rs Riemann solver
    \param hbc Hydrodynamic boundary conditions
   */
  ModularFluxCalculator
  (const SpatialReconstruction& sr,
   const RiemannSolver& rs,
   const HydroBoundaryConditions& hbc);

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
  const SpatialReconstruction& sr_;
  const RiemannSolver& rs_;
  const HydroBoundaryConditions& hbc_;
};

#endif // MODULAR_FLUX_CALCULATOR_HPP
