#ifndef DIFFERENT_BC_1D_HPP
#define DIFFERENT_BC_1D_HPP 1

#include "boundary_conditions_1d.hpp"

class DifferentBC: public BoundaryConditions1D
{
public:
  DifferentBC(const BoundaryConditions1D& left,
	      const BoundaryConditions1D& right);

  Extensive operator()
  (const SimulationState1D& ss,
   const EquationOfState& eos,
   const RiemannSolver& rs,
   const vector<double>& vertex_velocity,
   const bool side) const;
    
private:
    
  const BoundaryConditions1D& left_;
  const BoundaryConditions1D& right_;
};

#endif // DIFFERENT_BC_1D_HPP
