#ifndef LAGEU_HPP
#define LAGEU_HPP
#include "source/newtonian/two_dimensional/extensive_updater.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence_2.hpp"

class LagrangianExtensiveUpdater : public ExtensiveUpdater
{
public:
	LagrangianExtensiveUpdater(LagrangianFluxT const& fc);

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

}; 
#endif //LAGEU_HPP
