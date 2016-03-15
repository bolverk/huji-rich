/*! \file ColdFlowsExtensiveCalculator.hpp
\author Elad Steinberg
\brief Cold flows extensive variable updater, uses entropy for thermal energy
*/

#ifndef COLD_EXTENSIVE_UPDATER_HPP
#define COLD_EXTENSIVE_UPDATER_HPP 1

#include "ConditionExtensiveUpdater.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"
#include "interpolations/LinearGaussImproved.hpp"

//! \brief Cold flows extensive variable update
class ColdFlowsExtensiveCalculator : public ExtensiveUpdater
{
private:
	mutable ColdFlowsUpdate coldupdate_;
public:

	/*! \brief Class constructor
	  \param eos Equation of state
	  \param ghost The ghost point generator
	 */
	ColdFlowsExtensiveCalculator(EquationOfState const& eos, GhostPointGenerator const& ghost, LinearGaussImproved const& interp);

	void operator()
		(const vector<Extensive>& fluxes,
			const PhysicalGeometry& pg,
			const Tessellation& tess,
			const double dt,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			vector<Extensive>& extensives,
			double time,
			TracerStickerNames const& tracerstickernames) const;
};

#endif // COLD_EXTENSIVE_UPDATER_HPP
