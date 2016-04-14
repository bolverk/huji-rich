/*! \file flux_calculator_2d.hpp
  \author Almog Yalinewich
  \brief Base class for flux calculator
 */

#ifndef FLUX_CALCULATOR_HPP
#define FLUX_CALCULATOR_HPP 1

#include <vector>
#include "extensive.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell_2d.hpp"
#include "../common/equation_of_state.hpp"
#include "cache_data.hpp"

using std::vector;

//! \brief Base class for flux calculator
class FluxCalculator
{
public:

  /*! \brief Calculates fluxes
    \param tess Tessellation
    \param edge_velocities Velocities of the mesh generating points
    \param cells Computational cells
    \param extensives Extensive variables in each cell
    \param cd Cached data
    \param eos Equation of state
    \param time Time
    \param dt Time step
	\param tracerstickernames The names of the tracers and stickers
    \return List of fluxes on each edge
   */
  virtual vector<Extensive> operator()
  (const Tessellation& tess,
   const vector<Vector2D>& edge_velocities,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& extensives,
   const CacheData& cd,
   const EquationOfState& eos,
   const double time,
   const double dt,
	  TracerStickerNames const& tracerstickernames) const = 0;

  //! \brief Class destructor
  virtual ~FluxCalculator(void);
};

#endif // FLUX_CALCULATOR_HPP
