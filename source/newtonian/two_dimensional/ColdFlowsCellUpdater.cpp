#include "ColdFlowsCellUpdater.hpp"

ColdFlowsCellUpdate::ColdFlowsCellUpdate(const vector<pair<const SimpleCellUpdater::Condition*, 
	const SimpleCellUpdater::Action*> > sequence):scu_(SimpleCellUpdater(sequence)) {}

vector<ComputationalCell> ColdFlowsCellUpdate::operator()
(const Tessellation& tess,
	const PhysicalGeometry& pg,
	const EquationOfState& eos,
	vector<Extensive>& extensives,
	const vector<ComputationalCell>& old,
	const CacheData& cd) const
{
	vector<ComputationalCell> res = scu_(tess, pg, eos, extensives, old, cd);
	for (size_t i = 0; i < res.size(); ++i)
	{
		res[i].tracers["Entropy"] = eos.dp2s(res[i].density, res[i].pressure);
		extensives[i].tracers["Entropy"] = res[i].tracers["Entropy"] * cd.volumes[i] * res[i].density;
	}
	return res;
}
