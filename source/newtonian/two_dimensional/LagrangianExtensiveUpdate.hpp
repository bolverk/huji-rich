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
    \param eos Equation of state
    \param ghost Scheme to create ghost points
    \param interp Interpolation scheme
   */
	LagrangianExtensiveUpdate(LagrangianFlux const& lflux, EquationOfState const& eos, GhostPointGenerator
		const& ghost, LinearGaussImproved const& interp);

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
	EquationOfState const& eos_;
	GhostPointGenerator const& ghost_;
	LinearGaussImproved const& interp_;
};
#endif //LEUPDATE_HPP
