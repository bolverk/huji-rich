#ifndef LAGEU_HPP
#define LAGEU_HPP
#include "extensive_updater.hpp"
#include "condition_action_sequence_2.hpp"

class LagrangianExtensiveUpdater : public ExtensiveUpdater
{
public:
	LagrangianExtensiveUpdater(LagrangianFluxT const& fc,ExtensiveUpdater const& beu, EquationOfState const& eos);

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
	LagrangianFluxT const& fc_;
	ExtensiveUpdater const& beu_;
	EquationOfState const& eos_;
}; 
#endif //LAGEU_HPP
