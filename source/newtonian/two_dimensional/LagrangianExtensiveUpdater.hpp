#ifndef LAGEU_HPP
#define LAGEU_HPP
#include "extensive_updater.hpp"
#include "condition_action_sequence_2.hpp"

class LagrangianExtensiveUpdater : public ExtensiveUpdater
{
public:
	LagrangianExtensiveUpdater(LagrangianFlux const& fc,ExtensiveUpdater const& beu);

	void operator()
	(const vector<Extensive>& fluxes,
	const PhysicalGeometry& pg,
	const Tessellation& tess,
	const double dt,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	vector<Extensive>& extensives,
	double time, TracerStickerNames const& tracerstickersnames) const;
	
private:
	LagrangianFlux const& fc_;
	ExtensiveUpdater const& beu_;
}; 
#endif //LAGEU_HPP
