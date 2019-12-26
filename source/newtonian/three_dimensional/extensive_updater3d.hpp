/*! \file extensive_updater3d.hpp
\author Elad Steinberg
\brief Base class for extensive updater scheme
*/

#ifndef EXTENSIVE_UPDATER3D_HPP
#define EXTENSIVE_UPDATER3D_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "../common/equation_of_state.hpp"

using std::vector;

//! \brief Base class for extensive update scheme
class ExtensiveUpdater3D
{
public:

	/*! \brief Updates the extensive variables
	\param fluxes Fluxes
	\param tess Tessellation
	\param dt Time step
	\param cells Computational cells
	\param extensives Extensive variables
	\param tracerstickernames The names of the tracers and stickers
	\param time The time
	*/
	virtual void operator()(const vector<Conserved3D>& fluxes,const Tessellation3D& tess,
		const double dt,const vector<ComputationalCell3D>& cells,vector<Conserved3D>& extensives,double time,
		TracerStickerNames const& tracerstickernames, const vector<Vector3D>& edge_velocities,
		std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > const& interp_values) const = 0;

	//! \brief Class constructor
	virtual ~ExtensiveUpdater3D(void);
};

#endif // EXTENSIVE_UPDATER3D_HPP
