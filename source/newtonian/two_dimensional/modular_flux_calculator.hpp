#ifndef MODULAR_FLUX_CALCULATOR_HPP
#define MODULAR_FLUX_CALCULATOR_HPP 1

#include "flux_calculator_2d.hpp"
#include "spatial_reconstruction.hpp"
#include "../common/riemann_solver.hpp"
#include "GhostPointGenerator.hpp"

//! \brief Modular flux calculator
class ModularFluxCalculator: public FluxCalculator
{
public:

  /*! \brief Class constructor
    \param sr Interpolation
    \param rs Riemann solver
   */
  ModularFluxCalculator
  (const SpatialReconstruction& sr,
   const RiemannSolver& rs);

  vector<Extensive> operator()(const Tessellation& tess,const vector<Vector2D>& edge_velocities,
	  const vector<ComputationalCell>& cells,const vector<Extensive>& extensives,const CacheData& cd,
	  const EquationOfState& eos,const double time,const double dt,TracerStickerNames const& tracerstickernames) const;

private:
	const SpatialReconstruction& sr_;
	const RiemannSolver& rs_;
	mutable vector<pair<ComputationalCell, ComputationalCell> > interpolated_;
};

#endif // MODULAR_FLUX_CALCULATOR_HPP
