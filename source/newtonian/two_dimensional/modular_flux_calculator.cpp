#include "modular_flux_calculator.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"

ModularFluxCalculator::ModularFluxCalculator
(const SpatialReconstruction& sr,
 const RiemannSolver& rs,
 const HydroBoundaryConditions& hbc):
  sr_(sr), rs_(rs), hbc_(hbc),interpolated_(vector<pair<ComputationalCell, ComputationalCell> >()) {}

namespace
{
  pair<Vector2D, Vector2D> calc_parallel_normal(const Tessellation& tess,const Edge& edge)
  {
    const Vector2D p = Parallel(edge);
    if(edge.neighbors.first>=0 && edge.neighbors.first<tess.GetPointNo())
      return pair<Vector2D,Vector2D> (p,remove_parallel_component(edge.vertices.first - 
	  tess.GetMeshPoint(edge.neighbors.first),p));
    if(edge.neighbors.second>=0 && edge.neighbors.second<tess.GetPointNo())
      return pair<Vector2D,Vector2D> (p,remove_parallel_component(tess.GetMeshPoint(edge.neighbors.second) - 
	  edge.vertices.first,p));
    assert(false);				     
  }

  Extensive convert_conserved_to_extensive(const Conserved& conserved,const pair<ComputationalCell,ComputationalCell>& cells)
  {
    Extensive res;
    res.mass = conserved.Mass;
    res.momentum = conserved.Momentum;
    res.energy = conserved.Energy;
	res.tracers.reserve(cells.first.tracers.size());
	for (size_t i = 0; i < cells.first.tracers.size(); ++i)
	  res.tracers.insert
	    (pair<string, double>
	     ((cells.first.tracers.begin()+static_cast<int>(i))->first,
	      conserved.Mass*
	      (conserved.Mass>0 ? (cells.first.tracers.begin() + static_cast<int>(i))->second : 
	       (cells.second.tracers.begin() + static_cast<int>(i))->second)));
    return res;
  }

 }

vector<Extensive> ModularFluxCalculator::operator() (const Tessellation& tess,const vector<Vector2D>& edge_velocities,
	const vector<ComputationalCell>& cells,const vector<Extensive>& /*extensives*/,const CacheData& /*cd*/,
	const EquationOfState& eos,const double time,const double /*dt*/) const
{
	interpolated_.resize(tess.GetTotalSidesNumber(),
		pair<ComputationalCell, ComputationalCell>(cells[0], cells[0]));
	sr_(tess, cells, time,interpolated_);
  vector<bool> flags(static_cast<size_t>(tess.getAllEdges().size()), false);
  const vector<pair<size_t, Extensive> > boundary_conditions = hbc_(tess, cells);
  vector<Extensive> res(tess.getAllEdges().size());
  for(size_t i=0;i<boundary_conditions.size();++i)
  {
    const size_t index = boundary_conditions.at(i).first;
    flags.at(index) = true;
    res.at(index) = boundary_conditions.at(i).second;
  }
  for(size_t i=0;i<tess.getAllEdges().size();++i)
  {
    if(!flags.at(i))
	{
      flags.at(i) = true;
      const pair<Vector2D, Vector2D> p_n = calc_parallel_normal(tess,tess.getAllEdges().at(i));
      const Edge edge = tess.getAllEdges().at(i);
      const double speed = ScalarProd(p_n.second, edge_velocities.at(i)) / abs(p_n.second);
      const Primitive p_left = 
	convert_to_primitive
	(interpolated_.at(i).first, eos);
      const Primitive p_right = 
	convert_to_primitive
	(interpolated_.at(i).second, eos);
      res.at(i) = 
	convert_conserved_to_extensive
	(rotate_solve_rotate_back
	 (rs_,p_left,p_right,
	  speed,p_n.second,p_n.first),interpolated_.at(i));
    }
  }
  return res;
}
