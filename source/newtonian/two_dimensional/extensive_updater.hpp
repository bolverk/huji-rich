/*! \file extensive_updater.hpp
  \author Almog Yalinewich
  \brief Base class for extensive updater scheme
 */

#ifndef EXTENSIVE_UPDATER_HPP
#define EXTENSIVE_UPDATER_HPP 1

#include <vector>
#include "extensive.hpp"
#include "physical_geometry.hpp"
#include "../../tessellation/tessellation.hpp"
#include "cache_data.hpp"
#include "computational_cell_2d.hpp"

using std::vector;

//! \brief Base class for extensive update scheme
class ExtensiveUpdater
{
public:

	/*! \brief Updates the extensive variables
	  \param fluxes Fluxes
	  \param pg Physical geometry
	  \param tess Tessellation
	  \param dt Time step
	  \param cd Cache data
	  \param cells Computational cells
	  \param extensives Extensive variables
	  \param tracerstickernames The names of the tracers and stickers
	  \param time The time
	 */
	virtual void operator()
		(const vector<Extensive>& fluxes,
			const PhysicalGeometry& pg,
			const Tessellation& tess,
			const double dt,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			vector<Extensive>& extensives,
			double time,
			TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Class constructor
	virtual ~ExtensiveUpdater(void);
};

#endif // EXTENSIVE_UPDATER_HPP
