#include "condition_action_sequence_2.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"

ConditionActionSequence2::ConditionActionSequence2
(const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence,
	const vector<pair<const ConditionActionSequence::Condition*, const Action2*> >& sequence2,
	SpatialReconstruction const& interp):
	sequence_(sequence),sequence2_(sequence2),interp_(interp) {}

ConditionActionSequence2::~ConditionActionSequence2(void)
{
	/*
	for (size_t i = 0; i<sequence_.size(); ++i) {
		delete sequence_[i].first;
		delete sequence_[i].second;
	}
	for (size_t i = 0; i<sequence2_.size(); ++i) {
		delete sequence2_[i].first;
		delete sequence2_[i].second;
	}
	*/
}

namespace 
{
	Extensive choose_action
		(const Edge& edge,
			const Tessellation& tess,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos,
			const Vector2D& edge_velocity,
			const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence,
		const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*> >& sequence2,
			pair<ComputationalCell,ComputationalCell> const& edge_values)
	{
		for (size_t i = 0; i<sequence.size(); ++i) {
			const pair<bool, bool> flag_aux = (*sequence[i].first)
				(edge, tess, cells);
			if (flag_aux.first)
				return (*sequence[i].second)
				(edge, tess, edge_velocity, cells, eos, flag_aux.second);
		}
		for (size_t i = 0; i<sequence2.size(); ++i) {
			const pair<bool, bool> flag_aux = (*sequence2[i].first)
				(edge, tess, cells);
			if (flag_aux.first)
				return (*sequence2[i].second)
				(edge, tess, edge_velocity, cells, eos, flag_aux.second, edge_values);
		}
		throw "Error in ConditionActionSequence";
	}
}

vector<Extensive> ConditionActionSequence2::operator()
(const Tessellation& tess,
	const vector<Vector2D>& edge_velocities,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& /*extensives*/,
	const CacheData& /*cd*/,
	const EquationOfState& eos,
	const double time,
	const double /*dt*/) const
{
	vector<Extensive> res(tess.getAllEdges().size());
	vector<pair<ComputationalCell,ComputationalCell> > edge_values= interp_.operator()(tess, cells, time);
	for (size_t i = 0; i<tess.getAllEdges().size(); ++i)
		res[i] = choose_action
		(tess.getAllEdges()[i],
			tess,
			cells,
			eos,
			edge_velocities[i],
			sequence_,
			sequence2_,
			edge_values[i]);
	return res;
}


ConditionActionSequence2::Action2::~Action2(void) {}

RegularFlux2::RegularFlux2(const RiemannSolver& rs) :
	rs_(rs) {}

namespace {
	Extensive conserved_to_extensive
		(const Conserved& c, const ComputationalCell& cell)
	{
		Extensive res;
		res.mass = c.Mass;
		res.momentum = c.Momentum;
		res.energy = c.Energy;
		for (boost::container::flat_map<string, double>::const_iterator it =
			cell.tracers.begin();
			it != cell.tracers.end(); ++it)
			res.tracers[it->first] = (it->second)*c.Mass;
		return res;
	}
}

Extensive RegularFlux2::operator()
(const Edge& edge,
	const Tessellation& tess,
	const Vector2D& edge_velocity,
	const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,
	const bool /*aux*/,
	pair<ComputationalCell,ComputationalCell> const& edge_values) const
{
	const Vector2D p = normalize
		(edge.vertices.second -
			edge.vertices.first);
	const Vector2D n = normalize
		(tess.GetMeshPoint(edge.neighbors.second) -
			tess.GetMeshPoint(edge.neighbors.first));
	const double v =
		ScalarProd(n, edge_velocity);
	const Conserved c = rotate_solve_rotate_back
		(rs_,
			convert_to_primitive
			(edge_values.first, eos),
			convert_to_primitive
			(edge_values.second, eos),
			v, n, p);
	return conserved_to_extensive
		(c,
			c.Mass>0 ?
			edge_values.first :
			edge_values.second);
}

RigidWallFlux2::RigidWallFlux2
(const RiemannSolver& rs) :
	rs_(rs) {}

namespace {
	pair<Primitive, Primitive> rigid_wall_states
		(const Primitive& state,
			const Vector2D& p,
			const bool aux)
	{
		if (aux) {
			const Primitive left = state;
			const Primitive right = reflect(left, p);
			return pair<Primitive, Primitive>(left, right);
		}
		else {
			const Primitive right = state;
			const Primitive left = reflect(right, p);
			return pair<Primitive, Primitive>(left, right);
		}
	}
}

Extensive RigidWallFlux2::operator()
(const Edge& edge,
	const Tessellation& tess,
	const Vector2D& /*edge_velocity*/,
	const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,
	const bool aux,
	pair<ComputationalCell,ComputationalCell> const& edge_values) const
{
	if (aux)
		assert(edge.neighbors.first >= 0 && edge.neighbors.first<tess.GetPointNo());
	else
		assert(edge.neighbors.second >= 0 && edge.neighbors.second<tess.GetPointNo());
	const Vector2D p = normalize
		(edge.vertices.second - edge.vertices.first);
	const Vector2D n =
		normalize
		(remove_parallel_component
			(aux ?
				edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first) :
				tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
				p));
	const double v = 0;
	const pair<Primitive, Primitive> left_right =
		rigid_wall_states
		(convert_to_primitive
			(aux ? edge_values.first :edge_values.second,eos),
			p, aux);
	const Conserved c = rotate_solve_rotate_back
		(rs_,
			left_right.first,
			left_right.second,
			v, n, p);
	return conserved_to_extensive
		(c,
			aux ? edge_values.first : edge_values.second);
}

Ratchet::Ratchet(const RiemannSolver& rs,bool in) :	
  in_(in),
  wall_(RigidWallFlux2(rs)),
  free_(FreeFlowFlux(rs)) {}


Extensive Ratchet::operator()
(const Edge& edge,
	const Tessellation& tess,
	const Vector2D& edge_velocity,
	const vector<ComputationalCell>& cells,
	const EquationOfState& eos,
	const bool aux,
	const pair<ComputationalCell, ComputationalCell> & edge_values) const
{
	Vector2D n = aux ? tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first) :
		tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second);
	if (ScalarProd(n, cells[static_cast<size_t>(aux ? edge.neighbors.first : edge.neighbors.second)].velocity)*(2*static_cast<double>(in_)-1) < 0)
		return free_.operator()(edge,tess, edge_velocity, cells, eos, aux);
	else
		return wall_.operator()(edge, tess, edge_velocity, cells, eos, aux,edge_values);
}
