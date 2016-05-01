/*! \file SourceTerm.hpp
  \brief Abstract class for source terms
  \author Elad Steinberg
*/

#ifndef SOURCETERM_HPP
#define SOURCETERM_HPP 1
#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../../misc/utils.hpp"
#include "physical_geometry.hpp"
#include "extensive.hpp"
#include "computational_cell_2d.hpp"
#include "cache_data.hpp"

//! \brief Abstract class for external forces
class SourceTerm
{
public:
	/*!
	  \brief Calcualtes the change in conserved variables done on a cell from a source term
	  \param tess The tessellation
	  \param pg Physical geometry
	  \param cd Cache data
	  \param cells The hydrodynmic variables of the cell
	  \param fluxes Fluxes
	  \param point_velocities Velocities of the mesh generating points
	  \param t Time
	  \param tracerstickernames The names of the tracers and stickers
	  \return The flux of conserved variables
	*/
	virtual vector<Extensive> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& fluxes,
			const vector<Vector2D>& point_velocities,
			const double t,
			TracerStickerNames const& tracerstickernames) const = 0;

	virtual ~SourceTerm(void);
};

#endif //SOURCETERM_HPP
