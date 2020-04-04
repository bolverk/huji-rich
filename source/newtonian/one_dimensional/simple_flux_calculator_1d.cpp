#include "simple_flux_calculator_1d.hpp"

SimpleFluxCalculator1D::SimpleFluxCalculator1D
(const RiemannSolver& rs,
 const SpatialReconstruction1D& interp,
 const BoundaryConditions1D& bc):
  rs_(rs),
  interp_(interp),
  bc_(bc) {}

namespace{
  vector<Primitive> cc2primitives
    (const vector<ComputationalCell>& ccs,
     const EquationOfState& eos)
  {
    vector<Primitive> res(ccs.size());
    for(size_t i=0;i<res.size();++i){
      res.at(i).Density = ccs.at(i).density;
      res.at(i).Pressure = ccs.at(i).pressure;
      res.at(i).Velocity = ccs.at(i).velocity;
      res.at(i).Energy = eos.dp2e
	(ccs.at(i).density, ccs.at(i).pressure);
      res.at(i).SoundSpeed = eos.dp2c
	(ccs.at(i).density, ccs.at(i).pressure);
    }
    return res;
  }

  void copy_hydro_flux
    (const Conserved& hydro_flux,
     Extensive& flux)
  {
    flux.mass = hydro_flux.Mass;
    flux.momentum = hydro_flux.Momentum;
    flux.energy = hydro_flux.Energy;
  }

  void nullify_tracers(Extensive& flux)
  {
    for(size_t i=0;i<flux.tracers.size();++i)
      flux.tracers.at(i) = 0;
  }
}

vector<Extensive> SimpleFluxCalculator1D::operator()
(const SimulationState1D& ss,
 const vector<double>& vertex_velocity,
 const EquationOfState& eos,
 const double dt) const
{
  // Bulk
  vector<Extensive> res(ss.getVertices().size());
  for(size_t i=1;i<ss.getVertices().size()-1;++i){
    const Primitive left = interp_
      (ss.getVertices(),
       cc2primitives(ss.getCells(), eos),
       vertex_velocity.at(i),
       i,
       0,
       dt);
    const Primitive right = interp_
      (ss.getVertices(),
       cc2primitives(ss.getCells(), eos),
       vertex_velocity.at(i),
       i,
       1,
       dt);
    const Conserved hydro_flux = rs_(left, right, vertex_velocity.at(i));
    Extensive& flux = res.at(i);
    flux.mass = hydro_flux.Mass;
    flux.momentum = hydro_flux.Momentum;
    flux.energy = hydro_flux.Energy;

    const size_t source_index = hydro_flux.Mass>0 ? i-1 : i;
    for(size_t j=0;j<ss.getCells().at(i).tracers.size();++j)
      flux.tracers.push_back
	(res.at(i).mass*
	 ss.getCells().at(source_index).tracers.at(j));
  }

  // Boundary conditions
  const Conserved front = bc_(ss.getVertices(),
			      cc2primitives(ss.getCells(),eos),
			      rs_,
			      vertex_velocity,
			      0);
  copy_hydro_flux(front, res.front());
  res.front().tracers = ss.getCells().front().tracers;
  nullify_tracers(res.front());
  const Conserved back = bc_(ss.getVertices(),
			     cc2primitives(ss.getCells(),eos),
			     rs_,
			     vertex_velocity,
			     static_cast<int>(res.size()-1));
  copy_hydro_flux(back, res.back());
  res.back().tracers = ss.getCells().back().tracers;
  nullify_tracers(res.back());

  return res;			   
}
