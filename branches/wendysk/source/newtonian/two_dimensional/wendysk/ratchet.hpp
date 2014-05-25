#ifndef RATCHET_HPP
#define RATCHET_HPP 1

#include "../CustomEvolution.hpp"

class Ratchet: public CustomEvolution
{
public:

  enum DIRECTION {in, out};

  Ratchet(DIRECTION dir,
	  vector<Primitive> const& init_cells);

  bool flux_indifferent(void) const;

  Conserved CalcFlux(Tessellation const* tess,
		     vector<Primitive> const& cells,
		     double /*dt*/,
		     SpatialReconstruction* /*interp*/,
		     Edge const& edge,
		     Vector2D const& face_velocity,
		     RiemannSolver const& rs,
		     int index,
		     HydroBoundaryConditions const* hbc,
		     double /*time*/,
		     vector<vector<double> > const& /*tracers*/);

  Primitive UpdatePrimitive(vector<Conserved> const& /*intensives*/,
			    EquationOfState const* /*eos*/,
			    vector<Primitive>& /*cells*/,
			    int index);

private:
  const DIRECTION dir_;
  vector<Primitive> const& init_cells_;
};

#endif // RATCHET_HPP
