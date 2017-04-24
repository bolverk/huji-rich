#include "ConditionActionFlux1.hpp"

ConditionActionFlux1::ConditionActionFlux1(const vector<pair<const Condition3D*, const Action3D*> >& sequence,
	SpatialReconstruction3D const& interp) :
	sequence_(sequence),interp_(interp) {}

ConditionActionFlux1::~ConditionActionFlux1(void)
{}

namespace
{
	void choose_action(size_t face, const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const EquationOfState& eos, const Vector3D& face_velocity,
		const vector<pair<const ConditionActionFlux1::Condition3D*, const ConditionActionFlux1::Action3D*> >& sequence,
		Conserved3D &res, double time, TracerStickerNames const& tracerstickernames, std::pair<ComputationalCell3D, ComputationalCell3D>
		const& face_values)
	{
		for (size_t i = 0; i < sequence.size(); ++i)
		{
			const pair<bool, bool> flag_aux = (*sequence[i].first)
				(face, tess, cells, tracerstickernames);
			if (flag_aux.first)
				return (*sequence[i].second)
				(face, tess, face_velocity, cells, eos, flag_aux.second, res, time, tracerstickernames,face_values);
		}
		throw UniversalError("Error in ConditionActionFlux1");
	}
}

void ConditionActionFlux1::operator()(vector<Conserved3D> &fluxes, const Tessellation3D& tess, const vector<Vector3D>& face_velocities,
	const vector<ComputationalCell3D>& cells, const vector<Conserved3D>& extensives, const EquationOfState& eos,
	const double time, const double /*dt*/, TracerStickerNames const& tracerstickernames) const
{
	vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values;
	interp_(tess, cells, time, face_values, tracerstickernames);
	fluxes.resize(tess.GetTotalFacesNumber(), extensives[0]);
	size_t Nloop = fluxes.size();
	for (size_t i = 0; i < Nloop; ++i)
		choose_action(i, tess, cells, eos, face_velocities[i], sequence_, fluxes[i], time, tracerstickernames,face_values[i]);
}

ConditionActionFlux1::Condition3D::~Condition3D(void) {}

ConditionActionFlux1::Action3D::~Action3D(void) {}

RegularFlux3D::RegularFlux3D(const RiemannSolver3D& rs) :
	rs_(rs) {}

namespace
{
	void conserved_to_extensive
		(const Conserved3D& c, const ComputationalCell& cell, Conserved3D &res)
	{
		res.mass = c.mass;
		res.momentum = c.momentum;
		res.energy = c.energy;
		res.tracers.resize(cell.tracers.size());
		size_t N = cell.tracers.size();
		for (size_t i = 0; i < N; ++i)
			res.tracers[i] = cell.tracers[i] * c.mass;
	}
}

void RegularFlux3D::operator()(size_t face_index, const Tessellation3D& tess, const Vector3D& face_velocity,
	const vector<ComputationalCell3D>& /*cells*/, const EquationOfState& eos, const bool /*aux*/, Conserved3D &res,
	double /*time*/, TracerStickerNames const& tracerstickernames, std::pair<ComputationalCell3D, ComputationalCell3D>
	const& face_values) const
{
	const Vector3D normal = normalize(tess.Normal(face_index));
	RotateSolveBack3D(normal, face_values.first, face_values.second, face_velocity, rs_, res, eos, tracerstickernames);
}

RigidWallFlux3D::RigidWallFlux3D(const RiemannSolver3D& rs) : rs_(rs) {}

namespace
{
	void reflect(ComputationalCell3D & cell, Vector3D const& normal)
	{
		cell.velocity -= 2 * ScalarProd(cell.velocity, normal)*normal;
	}

	void rigid_wall_states(std::pair<ComputationalCell3D, ComputationalCell3D> &res,
		Vector3D const& normal, const bool aux)
	{
		if (aux)
		{
			res.second = res.first;
			reflect(res.second, normal);
		}
		else
		{
			res.first = res.second;
			reflect(res.first, normal);
		}
	}
}

void RigidWallFlux3D::operator()(size_t face_index, const Tessellation3D& tess, const Vector3D& face_velocity,
	const vector<ComputationalCell3D>& /*cells*/, const EquationOfState& eos, const bool aux, Conserved3D &res,
	double /*time*/, TracerStickerNames const& tracerstickernames, std::pair<ComputationalCell3D, ComputationalCell3D>
	const& face_values) const
{
	const Vector3D normal = normalize(tess.Normal(face_index));
	std::pair<ComputationalCell3D, ComputationalCell3D>	rigid_states(face_values);
	rigid_wall_states(rigid_states, normal, aux);
	RotateSolveBack3D(normal, rigid_states.first, rigid_states.second, face_velocity, rs_, res, eos, tracerstickernames);
}

FreeFlowFlux3D::FreeFlowFlux3D(const RiemannSolver3D& rs) : rs_(rs) {}

void FreeFlowFlux3D::operator()(size_t face_index, const Tessellation3D& tess, const Vector3D& face_velocity,
	const vector<ComputationalCell3D>& /*cells*/, const EquationOfState& eos, const bool aux, Conserved3D &res,
	double /*time*/, TracerStickerNames const& tracerstickernames, std::pair<ComputationalCell3D, ComputationalCell3D>
	const& face_values) const
{
	const Vector3D normal = normalize(tess.Normal(face_index));
	std::pair<ComputationalCell3D, ComputationalCell3D>	states(face_values);
	if (aux)
		states.second = states.first;
	else
		states.first = states.second;
	RotateSolveBack3D(normal, states.first, states.second, face_velocity, rs_, res, eos, tracerstickernames);
}

IsBoundaryFace3D::IsBoundaryFace3D(void) {}

pair<bool, bool> IsBoundaryFace3D::operator()(size_t face_index, const Tessellation3D& tess,
	const vector<ComputationalCell3D>& /*cells*/, TracerStickerNames const& /*tracerstickernames*/) const
{
	if (!tess.BoundaryFace(face_index))
		return pair<bool, bool>(false, false);
	if (tess.GetFaceNeighbors(face_index).first > tess.GetPointNo())
		return pair<bool, bool>(true, false);
	else
		return pair<bool, bool>(true, true);
}

IsBulkFace3D::IsBulkFace3D(void) {}

pair<bool, bool> IsBulkFace3D::operator()(size_t face_index, const Tessellation3D& tess,
	const vector<ComputationalCell3D>& /*cells*/, TracerStickerNames const& /*tracerstickernames*/)const
{
	if (tess.BoundaryFace(face_index))
		return pair<bool, bool>(false, false);
	else
		return pair<bool, bool>(true, false);
}

RegularSpecialEdge3D::RegularSpecialEdge3D(const string& sticker_name) :
	sticker_name_(sticker_name) {}

pair<bool, bool> RegularSpecialEdge3D::operator()(size_t face_index, const Tessellation3D& tess,
	const vector<ComputationalCell3D>& cells, TracerStickerNames const& tracerstickernames)const
{
	if (safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).first).stickers, tracerstickernames.sticker_names,
		sticker_name_))
	{
		if (safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).second).stickers, tracerstickernames.sticker_names,
			sticker_name_))
			return pair<bool, bool>(false, false);
		return pair<bool, bool>(true, false);
	}
	if (safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).second).stickers, tracerstickernames.sticker_names,
		sticker_name_))
		return pair<bool, bool>(true, true);
	return pair<bool, bool>(false, false);
}
