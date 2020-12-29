/*! \file cell_updater_2d.hpp
  \author Almog Yalinewich
  \brief Base class for cell update scheme
 */

#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include <vector>
#include "computational_cell_2d.hpp"
#include "../../tessellation/tessellation.hpp"
#include "physical_geometry.hpp"
#include "../common/equation_of_state.hpp"
#include "extensive.hpp"
#include "cache_data.hpp"

using std::vector;

//! \brief Base class for cell update scheme
class CellUpdater
{
public:

	/*! \brief Calculates the computational cells
	  \param tess Tessellation
	  \param pg Physical geometry
	  \param eos Equation of state
	  \param extensives Extensive variables
	  \param old Old computational cells
	  \param cd Cached data
	  \param tracerstickername The names of the tracers and stickers
	  \param time Time
	  \return List of computational cells
	 */
	virtual vector<ComputationalCell> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			vector<Extensive>& extensives,
			const vector<ComputationalCell>& old,
			const CacheData& cd,
			const TracerStickerNames& tracerstickername,
			double time) const = 0;

	//! \brief Class destructor
	virtual ~CellUpdater(void);
};

/*!
\brief Calculates velocity from extensive in SR
\param cell The extensive variable
\param G The adiabatic index
\return The velocity
*/
double GetVelocity(Extensive const& cell, double G);

#endif // CELL_UPDATER_HPP
