/*! \brief different_bc.hpp
  \author Almog Yalinewich
  \brief A class to combine two different boundary conditions
 */

#ifndef DIFFERENT_BC_1D_HPP
#define DIFFERENT_BC_1D_HPP 1

#include "boundary_conditions_1d.hpp"

//! \brief A combination of two different boundary conditions
class DifferentBC: public BoundaryConditions1D
{
public:

  /*! \brief Class constructor
    \param left Left boundary condition
    \param right Right boundary condtion
   */
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
