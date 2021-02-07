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
#ifdef RICH_DEBUG
			try
			{
#endif
				if (flag_aux.first)
					return (*sequence[i].second)
					(face, tess, face_velocity, cells, eos, flag_aux.second, res, time, tracerstickernames, face_values);
#ifdef RICH_DEBUG
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Face number", face);
				size_t n0 = tess.GetFaceNeighbors(face).first;
				size_t n1 = tess.GetFaceNeighbors(face).second;
				eo.AddEntry("Left cell density", cells[n0].density);
				eo.AddEntry("Left cell pressure", cells[n0].pressure);
				eo.AddEntry("Left cell Vx", cells[n0].velocity.x);
				eo.AddEntry("Left cell Vy", cells[n0].velocity.y);
				eo.AddEntry("Left cell Vz", cells[n0].velocity.z);
				eo.AddEntry("Left cell internal energy", cells[n0].internal_energy);
				eo.AddEntry("Left cell ID", cells[n0].ID);
				eo.AddEntry("Right cell density", cells[n1].density);
				eo.AddEntry("Right cell pressure", cells[n1].pressure);
				eo.AddEntry("Right cell Vx", cells[n1].velocity.x);
				eo.AddEntry("Right cell Vy", cells[n1].velocity.y);
				eo.AddEntry("Right cell Vz", cells[n1].velocity.z);
				eo.AddEntry("Right cell internal energy", cells[n1].internal_energy);
				eo.AddEntry("Right cell ID", cells[n1].ID);
				throw eo;
			}
#endif
		}
		throw UniversalError("Error in ConditionActionFlux1");
	}
}

std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > ConditionActionFlux1::operator()(vector<Conserved3D> &fluxes, const Tessellation3D& tess, const vector<Vector3D>& face_velocities,
												   const vector<ComputationalCell3D>& cells, const vector<Conserved3D>& /*extensives*/, const EquationOfState& eos,
	const double time, const double /*dt*/, TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < sequence_.size(); ++i)
		sequence_[i].second->Reset();
	vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values;
	interp_(tess, cells, time, face_values, tracerstickernames);
	fluxes.resize(tess.GetTotalFacesNumber());
	size_t Nloop = fluxes.size();
	for (size_t i = 0; i < Nloop; ++i)
	{
		if (face_values[i].first.density < 0 ||  face_values[i].first.internal_energy < 0
			|| face_values[i].second.density < 0  || face_values[i].second.internal_energy < 0)
		{
			UniversalError eo("Bad input to flux calculator");
			eo.AddEntry("Face", static_cast<double>(i));
			eo.AddEntry("Face neigh 0", static_cast<double>(tess.GetFaceNeighbors(i).first));
			eo.AddEntry("Face neigh 1", static_cast<double>(tess.GetFaceNeighbors(i).second));
			eo.AddEntry("First input Density", face_values[i].first.density);
			eo.AddEntry("First input pressure", face_values[i].first.pressure);
			eo.AddEntry("First input internal energy", face_values[i].first.internal_energy);
			eo.AddEntry("First input vx", face_values[i].first.velocity.x);
			eo.AddEntry("First input vy", face_values[i].first.velocity.y);
			eo.AddEntry("First input vz", face_values[i].first.velocity.z);
			eo.AddEntry("Second input Density", face_values[i].second.density);
			eo.AddEntry("Second input pressure", face_values[i].second.pressure);
			eo.AddEntry("Second input internal energy", face_values[i].second.internal_energy);
			eo.AddEntry("Second input vx", face_values[i].second.velocity.x);
			eo.AddEntry("Second input vy", face_values[i].second.velocity.y);
			eo.AddEntry("Second input vz", face_values[i].second.velocity.z);
			eo.AddEntry("First cell Density", cells[tess.GetFaceNeighbors(i).first].density);
			eo.AddEntry("First cell pressure", cells[tess.GetFaceNeighbors(i).first].pressure);
			eo.AddEntry("First cell internal energy", cells[tess.GetFaceNeighbors(i).first].internal_energy);
			eo.AddEntry("First cell vx", cells[tess.GetFaceNeighbors(i).first].velocity.x);
			eo.AddEntry("First cell vy", cells[tess.GetFaceNeighbors(i).first].velocity.y);
			eo.AddEntry("First cell vz", cells[tess.GetFaceNeighbors(i).first].velocity.z);
			eo.AddEntry("Second cell Density", cells[tess.GetFaceNeighbors(i).second].density);
			eo.AddEntry("Second cell pressure", cells[tess.GetFaceNeighbors(i).second].pressure);
			eo.AddEntry("Second cell internal energy", cells[tess.GetFaceNeighbors(i).second].internal_energy);
			eo.AddEntry("Second cell vx", cells[tess.GetFaceNeighbors(i).second].velocity.x);
			eo.AddEntry("Second cell vy", cells[tess.GetFaceNeighbors(i).second].velocity.y);
			eo.AddEntry("Second cell vz", cells[tess.GetFaceNeighbors(i).second].velocity.z);
			eo.AddEntry("Face vx", face_velocities[i].x);
			eo.AddEntry("Face vy", face_velocities[i].y);
			eo.AddEntry("Face vz", face_velocities[i].z);
			throw eo;
		}
		choose_action(i, tess, cells, eos, face_velocities[i], sequence_, fluxes[i], time, tracerstickernames, face_values[i]);
	}
	return face_values;
}

ConditionActionFlux1::Condition3D::~Condition3D(void) {}

ConditionActionFlux1::Action3D::~Action3D(void) {}

RegularFlux3D::RegularFlux3D(const RiemannSolver3D& rs) :
	rs_(rs) {}

/*namespace
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
}*/

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
	if (*safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).first).stickers.begin(), tracerstickernames.sticker_names.begin(),
		tracerstickernames.sticker_names.end(), sticker_name_))
	{
		if (*safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).second).stickers.begin(), tracerstickernames.sticker_names.begin(),
			tracerstickernames.sticker_names.end(), sticker_name_))
			return pair<bool, bool>(false, false);
		return pair<bool, bool>(true, false);
	}
	if (*safe_retrieve(cells.at(tess.GetFaceNeighbors(face_index).second).stickers.begin(), tracerstickernames.sticker_names.begin(),
		tracerstickernames.sticker_names.end(), sticker_name_))
		return pair<bool, bool>(true, true);
	return pair<bool, bool>(false, false);
}

LagrangianFlux3D::LagrangianFlux3D(const LagrangianHLLC3D & rs, const LagrangianHLLC3D & rs2, 
	LagrangianCriteria3D const & criteria):ws_(std::vector<double>()), edge_vel_(std::vector<double>()),
	Lag_calc_(std::vector<bool>()), rs_(rs), rs2_(rs2), criteria_(criteria) {}

void LagrangianFlux3D::operator()(size_t face_index, const Tessellation3D & tess, const Vector3D & face_velocity, 
	const vector<ComputationalCell3D>& cells, const EquationOfState & eos, const bool aux, Conserved3D & res, double time, 
	TracerStickerNames const & tracerstickernames, std::pair<ComputationalCell3D, ComputationalCell3D> const & face_values) const
{
	size_t N = tess.GetTotalFacesNumber();
	ws_.resize(N, 0.0);
	edge_vel_.resize(N, 0.0);
	Lag_calc_.resize(N, false);
	const Vector3D normal = normalize(tess.Normal(face_index));
	if (criteria_(face_index, tess, face_velocity, cells, eos, aux, face_values, time, tracerstickernames))
	{
		RotateSolveBack3D(normal, face_values.first, face_values.second, face_velocity, rs_, res, eos, tracerstickernames);
		ws_[face_index] = rs_.ws;
		Lag_calc_[face_index] = true;
	}
	else
	{
		RotateSolveBack3D(normal, face_values.first, face_values.second, face_velocity, rs2_, res, eos, tracerstickernames);
		ws_[face_index] = 0;
		Lag_calc_[face_index] = false;
	}
	edge_vel_[face_index] = ScalarProd(normal, face_velocity);
}

void LagrangianFlux3D::Reset(void) const
{
	ws_.assign(ws_.size(), 0);
	edge_vel_.assign(edge_vel_.size(), 0);
	Lag_calc_.assign(Lag_calc_.size(), false);
}

LagrangianFlux3D::LagrangianCriteria3D::~LagrangianCriteria3D(){}
