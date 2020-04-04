#include "different_bc.hpp"

DifferentBC::DifferentBC
(const BoundaryConditions1D& left,
 const BoundaryConditions1D& right):
  left_(left), right_(right) {}

Extensive DifferentBC::operator()
(const SimulationState1D& ss,
 const EquationOfState& eos,
 const RiemannSolver& rs,
 const vector<double>& vertex_velocity,
 const bool side) const
{
  const BoundaryConditions1D& bc =
    side ? right_ : left_;
  return bc(ss,eos,rs,vertex_velocity,side);
}
