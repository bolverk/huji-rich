#include "default_cell_updater.hpp"
#include "../../misc/utils.hpp"

DefaultCellUpdater::DefaultCellUpdater(void):entropy_index_(9999999){}

void DefaultCellUpdater::operator()(vector<ComputationalCell3D> &res, EquationOfState const& eos,
	const Tessellation3D& tess, vector<Conserved3D>& extensives, 
	TracerStickerNames const& tracerstickernames) const
{
	entropy_index_ = tracerstickernames.tracer_names.size();
	vector<string>::const_iterator it = binary_find(tracerstickernames.tracer_names.begin(),
		tracerstickernames.tracer_names.end(), string("Entropy"));
	if (it != tracerstickernames.tracer_names.end())
		entropy_index_ = static_cast<size_t>(it - tracerstickernames.tracer_names.begin());
	size_t Nloop = tess.GetPointNo();
	size_t Ntracers = res.at(0).tracers.size();
	res.resize(Nloop);
	for (size_t i = 0; i < Nloop; ++i)
	{
		Conserved3D& extensive = extensives[i];
		const double vol = tess.GetVolume(i);
		res[i].density = extensive.mass / vol;
		res[i].velocity = extensive.momentum / extensive.mass;
		double energy = extensive.internal_energy / extensive.mass;
		if (energy < 0)
		{
			energy = extensive.energy / extensive.mass - 0.5*ScalarProd(res[i].velocity, res[i].velocity);
			extensive.internal_energy = energy*extensive.mass;
		}
		extensive.energy = extensive.mass*(energy + 0.5*ScalarProd(res[i].velocity, res[i].velocity));
		//const double energy = extensive.energy / extensive.mass - 0.5*ScalarProd(res[i].velocity,res[i].velocity);
		res[i].stickers = res[i].stickers;
		for (size_t j = 0; j < Ntracers; ++j)
			res[i].tracers[j] = extensive.tracers[j] / extensive.mass;
		if (res[i].density < 0 || energy < 0)
		{
			UniversalError eo("Negative quantity in cell update");
			eo.AddEntry("Cell index", static_cast<double>(i));
			eo.AddEntry("Cell mass", extensive.mass);
			eo.AddEntry("Cell x momentum", extensive.momentum.x);
			eo.AddEntry("Cell y momentum", extensive.momentum.y);
			eo.AddEntry("Cell z momentum", extensive.momentum.z);
			eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
			eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
			eo.AddEntry("Cell z location", tess.GetMeshPoint(i).z);
			eo.AddEntry("Cell volume", vol);
			eo.AddEntry("Cell energy", extensive.energy);
			eo.AddEntry("Cell thermal energy per unit mass", energy);
			throw eo;
		}
		res[i].pressure = eos.de2p(res[i].density, energy, res[i].tracers, tracerstickernames.tracer_names);
		res[i].internal_energy = energy;
		if (entropy_index_ < res[i].tracers.size())
		{
			res[i].tracers[entropy_index_] = eos.dp2s(res[i].density, res[i].pressure, res[i].tracers,
				tracerstickernames.tracer_names);
			extensive.tracers[entropy_index_] = res[i].tracers[entropy_index_] * extensive.mass;
		}
	}
}
