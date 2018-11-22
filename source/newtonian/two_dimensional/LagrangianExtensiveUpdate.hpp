#ifndef LEUPDATE_HPP
#define LEUPDATE_HPP 1
#include "extensive_updater.hpp"
#include "condition_action_sequence_2.hpp"
#include "ConditionExtensiveUpdater.hpp"
#include "interpolations/LinearGaussImproved.hpp"

/*! \brief Updates extensive in such a way that minimises advection
 */
class LagrangianExtensiveUpdate : public ExtensiveUpdater
{
public:
  /*! \brief Class constructor
    \param lflux Flux calculator
    \param ghost Scheme to create ghost points
   */
	LagrangianExtensiveUpdate(LagrangianFlux const& lflux, GhostPointGenerator
		const& ghost);

	void operator()
		(const vector<Extensive>& fluxes,
			const PhysicalGeometry& pg,
			const Tessellation& tess,
			const double dt,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			vector<Extensive>& extensives,
			double time, TracerStickerNames const& tracerstickernames) const;
private:
	LagrangianFlux const& lflux_;
	GhostPointGenerator const& ghost_;
};
#endif //LEUPDATE_HPP
