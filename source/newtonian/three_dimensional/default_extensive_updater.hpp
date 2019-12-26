/*! \file default_extensive_updater.hpp
  \brief Generates extensive conserved variables
  \author Almog Yalinewich
 */

#ifndef DEFAULT_EXTENSIVE_UPDATER_HPP
#define DEFAULT_EXTENSIVE_UPDATER_HPP 1

#include "conserved_3d.hpp"
#include "extensive_updater3d.hpp"

//! \brief Generates a list of conserved variables
class DefaultExtensiveUpdater: public ExtensiveUpdater3D
{
public:

  /*! \brief Class constructor
   */
	DefaultExtensiveUpdater(void);

	void operator()(const vector<Conserved3D>& fluxes, const Tessellation3D& tess,
		const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double time,
		TracerStickerNames const& tracerstickernames, const vector<Vector3D>& edge_velocities,
		std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > const& interp_values) const;
};

#endif // DEFAULT_EXTENSIVE_UPDATER_HPP
