#include "condition_action_sequence_2.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"

ConditionActionSequence2::ConditionActionSequence2
(const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence,
	const vector<pair<const ConditionActionSequence::Condition*, const Action2*> >& sequence2,
	SpatialReconstruction const& interp):
	sequence_(sequence),sequence2_(sequence2),interp_(interp),edge_values_(vector<pair<ComputationalCell,
	ComputationalCell> >()){}

ConditionActionSequence2::~ConditionActionSequence2(void)
{}

namespace 
{
	pair<Vector2D, Vector2D> calc_parallel_normal(const Tessellation& tess, const Edge& edge)
	{
		const Vector2D p = Parallel(edge);
		if (edge.neighbors.first >= 0 && edge.neighbors.first < tess.GetPointNo())
			return pair<Vector2D, Vector2D>(p, remove_parallel_component(edge.vertices.first -
				tess.GetMeshPoint(edge.neighbors.first), p));
		if (edge.neighbors.second >= 0 && edge.neighbors.second < tess.GetPointNo())
			return pair<Vector2D, Vector2D>(p, remove_parallel_component(tess.GetMeshPoint(edge.neighbors.second) -
				edge.vertices.first, p));
		assert(false);
	}
	
	Extensive convert_conserved_to_extensive(const Conserved& conserved, const pair<ComputationalCell, ComputationalCell>& cells)
	{
		Extensive res;
		res.mass = conserved.Mass;
		res.momentum = conserved.Momentum;
		res.energy = conserved.Energy;
		res.tracers.resize(cells.first.tracers.size());
		size_t N = res.tracers.size();
		for (size_t i = 0; i < N; ++i)
			res.tracers[i] = conserved.Mass*(conserved.Mass>0 ? cells.first.tracers[i] : cells.second.tracers[i]);
		return res;
	}
	
	void choose_action
		(const Edge& edge,
			const Tessellation& tess,
			const vector<ComputationalCell>& cells,
			const EquationOfState& eos,
			const Vector2D& edge_velocity,
			const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >& sequence,
		const vector<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence2::Action2*> >& sequence2,
			pair<ComputationalCell,ComputationalCell> const& edge_values,Extensive &res,double time,
			TracerStickerNames const& tracerstickernames,size_t index)
	{
		for (size_t i = 0; i<sequence.size(); ++i) 
		{
			const pair<bool, bool> flag_aux = (*sequence[i].first)
				(edge, tess, cells,tracerstickernames);
			if (flag_aux.first)
			{
				(*sequence[i].second)(edge, tess, edge_velocity, cells, eos, flag_aux.second,res,time,tracerstickernames);
				return;
			}
		}
		for (size_t i = 0; i<sequence2.size(); ++i) 
		{
			const pair<bool, bool> flag_aux = (*sequence2[i].first)
				(edge, tess, cells,tracerstickernames);
			if (flag_aux.first)
			{
				(*sequence2[i].second)(edge,index, tess, edge_velocity, cells, eos, flag_aux.second, edge_values,res,time,
					tracerstickernames);
				return;
			}
		}
		throw UniversalError("Error in ConditionActionSequence");
	}
}

vector<Extensive> ConditionActionSequence2::operator()
(const Tessellation& tess,
	const vector<Vector2D>& edge_velocities,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& extensives,
	const CacheData& cd,
	const EquationOfState& eos,
	const double time,
	const double /*dt*/,
	TracerStickerNames const& tracerstickernames) const
{
	edge_values_.resize(static_cast<size_t>(tess.GetTotalSidesNumber()),
		pair<ComputationalCell, ComputationalCell>(cells[0], cells[0]));
	interp_.operator()(tess, cells, time,edge_values_,tracerstickernames,cd);
	vector<Extensive> res(tess.getAllEdges().size(), extensives[0]);
	for (size_t i = 0; i<tess.getAllEdges().size(); ++i)
		choose_action
		(tess.getAllEdges()[i],
			tess,
			cells,
			eos,
			edge_velocities[i],
			sequence_,
			sequence2_,
			edge_values_[i],res[i],time,tracerstickernames,i);
	return res;
}


ConditionActionSequence2::Action2::~Action2(void) {}

RegularFlux2::RegularFlux2(const RiemannSolver& rs) :
	rs_(rs) {}

namespace 
{
	void conserved_to_extensive
		(const Conserved& c, const ComputationalCell& cell,Extensive &res)
	{
		res.mass = c.Mass;
		res.momentum = c.Momentum;
		res.energy = c.Energy;
		res.tracers.resize(cell.tracers.size());
		size_t N = cell.tracers.size();
		for (size_t i = 0; i < N; ++i)
			res.tracers[i] = cell.tracers[i] * c.Mass;
	}
}

void RegularFlux2::operator()
(const Edge& edge,
	const size_t /*index*/,
	const Tessellation& tess,
	const Vector2D& edge_velocity,
	const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,
	const bool /*aux*/,
	pair<ComputationalCell,ComputationalCell> const& edge_values,
	Extensive &res,double /*time*/,
	TracerStickerNames const& tracerstickernames) const
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
			(edge_values.first, eos,tracerstickernames),
			convert_to_primitive
			(edge_values.second, eos,tracerstickernames),
			v, n, p);
	conserved_to_extensive(c,c.Mass>0 ?	edge_values.first :	edge_values.second,res);
}

RigidWallFlux2::RigidWallFlux2
(const RiemannSolver& rs) :
	rs_(rs) {}

namespace 
{
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

void RigidWallFlux2::operator()
(const Edge& edge,
	const size_t /*index*/,
	const Tessellation& tess,
	const Vector2D& /*edge_velocity*/,
	const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,
	const bool aux,
	pair<ComputationalCell,ComputationalCell> const& edge_values,
	Extensive &res,double /*time*/,
	TracerStickerNames const& tracerstickernames) const
{
#ifndef RICH_MPI
	if (aux)
		assert(edge.neighbors.first >= 0 && edge.neighbors.first<tess.GetPointNo());
	else
		assert(edge.neighbors.second >= 0 && edge.neighbors.second<tess.GetPointNo());
#endif
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
			(aux ? edge_values.first :edge_values.second,eos,tracerstickernames),
			p, aux);
	const Conserved c = rotate_solve_rotate_back
		(rs_,
			left_right.first,
			left_right.second,
			v, n, p);
	conserved_to_extensive(c,aux ? edge_values.first : edge_values.second,res);
}

Ratchet::Ratchet(const RiemannSolver& rs,bool in) :	
  in_(in),
  wall_(RigidWallFlux2(rs)),
  free_(FreeFlowFlux(rs)) {}


void Ratchet::operator()
(const Edge& edge,
	const size_t index,
	const Tessellation& tess,
	const Vector2D& edge_velocity,
	const vector<ComputationalCell>& cells,
	const EquationOfState& eos,
	const bool aux,
	const pair<ComputationalCell, ComputationalCell> & edge_values,
	Extensive &res,double time,
	TracerStickerNames const& tracerstickernames) const
{
	Vector2D n = aux ? tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first) :
		tess.GetMeshPoint(edge.neighbors.first) - tess.GetMeshPoint(edge.neighbors.second);
	if (ScalarProd(n, cells[static_cast<size_t>(aux ? edge.neighbors.first : edge.neighbors.second)].velocity)*(2*static_cast<double>(in_)-1) < 0)
		free_.operator()(edge,tess, edge_velocity, cells, eos, aux,res,time,tracerstickernames);
	else
		wall_.operator()(edge,index, tess, edge_velocity, cells, eos, aux,edge_values,res,time,tracerstickernames);
}

LagrangianFlux::LagrangianFlux(const LagrangianHLLC& rs):rs_(rs),ws_(vector<double>()),
	edge_vel_(vector<double>()){}
	
void LagrangianFlux::operator()(const Edge& edge,const size_t index,const Tessellation& tess,const Vector2D& edge_velocity,const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,const bool /*aux*/,const pair<ComputationalCell, ComputationalCell> & edge_values,Extensive &res, double /*time*/,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = static_cast<size_t>(tess.GetTotalSidesNumber());
	ws_.resize(N,0.0);
	edge_vel_.resize(N,0.0);
	const pair<Vector2D, Vector2D> p_n = calc_parallel_normal(tess, edge);
	const double speed = ScalarProd(p_n.second, edge_velocity) / abs(p_n.second);
	const Primitive p_left =convert_to_primitive(edge_values.first, eos, tracerstickernames);
	const Primitive p_right =convert_to_primitive(edge_values.second, eos, tracerstickernames);
	res =	convert_conserved_to_extensive(rotate_solve_rotate_back(rs_, p_left, p_right,speed, p_n.second, p_n.first), edge_values);
	ws_[index] = rs_.energy;
	edge_vel_[index] = speed;
}

LagrangianFluxT::LagrangianFluxT(const LagrangianHLLC& rs,const LagrangianHLLC& rs2):rs_(rs),rs2_(rs2),ws_(vector<double>()),
	edge_vel_(vector<double>()){}
	
void LagrangianFluxT::operator()(const Edge& edge,const size_t index,const Tessellation& tess,const Vector2D& edge_velocity,const vector<ComputationalCell>& /*cells*/,
	const EquationOfState& eos,const bool /*aux*/,const pair<ComputationalCell, ComputationalCell> & edge_values,Extensive &res, double /*time*/,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = static_cast<size_t>(tess.GetTotalSidesNumber());
	ws_.resize(N,0.0);
	edge_vel_.resize(N,0.0);
	const pair<Vector2D, Vector2D> p_n = calc_parallel_normal(tess, edge);
	const double speed = ScalarProd(p_n.second, edge_velocity) / abs(p_n.second);
	const Primitive p_left =convert_to_primitive(edge_values.first, eos, tracerstickernames);
	const Primitive p_right =convert_to_primitive(edge_values.second, eos, tracerstickernames);
	//if(edge_values.first.tracers.size()>14&&edge_values.first.tracers[14]<2.5e8&&edge_values.second.tracers[14]<2.5e8)
	if (edge.neighbors.first > tess.GetPointNo() || edge.neighbors.second > tess.GetPointNo())
	{
		res =	convert_conserved_to_extensive(rotate_solve_rotate_back(rs2_, p_left, p_right,speed, p_n.second, p_n.first), edge_values);
		ws_[index] = 0;
	}
	else
	{
		res =	convert_conserved_to_extensive(rotate_solve_rotate_back(rs_, p_left, p_right,speed, p_n.second, p_n.first), edge_values);
		ws_[index] = rs_.energy;
	}
	edge_vel_[index] = speed;
}