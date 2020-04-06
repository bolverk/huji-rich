#include "simple_flux_calculator_1d.hpp"
#include "flux_conversion.hpp"

SimpleFluxCalculator1D::SimpleFluxCalculator1D
(const RiemannSolver& rs,
 const SpatialReconstruction1D& interp,
 const BoundaryConditions1D& bc):
  rs_(rs),
  interp_(interp),
  bc_(bc) {}

vector<Extensive> SimpleFluxCalculator1D::operator()
(const SimulationState1D& ss,
 const vector<double>& vertex_velocity,
 const EquationOfState& eos,
 const double dt) const
{
  // Bulk
  vector<Extensive> res(ss.getVertices().size());
  const vector<Primitive> primitives = ccs2primitives
    (ss.getCells(), eos);
  for(size_t i=1;i<ss.getVertices().size()-1;++i){
    const Primitive left = interp_
      (ss.getVertices(),
       primitives,
       vertex_velocity.at(i),
       i,
       0,
       dt);
    const Primitive right = interp_
      (ss.getVertices(),
       primitives,
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
  res.front() = bc_(ss,
		    eos,
		    rs_,
		    vertex_velocity,
		    false);
  res.back() = bc_(ss,
		   eos,
		   rs_,
		   vertex_velocity,
		   true);
  return res;			   
}
