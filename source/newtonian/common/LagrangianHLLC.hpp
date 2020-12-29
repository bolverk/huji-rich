/*! \file LagrangianHLLC.hpp
\brief HLLC riemann solver that keeps track of transfered internal energy
\author Elad Steinberg
*/

#ifndef LAGRANGIANHLLC_HPP
#define LAGRANGIANHLLC_HPP 1

#include "riemann_solver.hpp"
#include "equation_of_state.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"

//! \brief LagrangianHLLC Riemann solver for an Eulerian grid
class LagrangianHLLC : public RiemannSolver
{
private:
	EquationOfState const& eos_;
	const bool massflux_,iter_;
public:
  /*! \brief Class constructor
    \param massflux Whether to apply correction for mass flux
    \param eos Equation of state
    \param iter Perform multiple iterations
   */
  explicit LagrangianHLLC(EquationOfState const& eos, bool massflux=true,bool iter=false);

	Conserved operator()
		(Primitive const& left,
			Primitive const& right,
			double velocity) const;

  //! \brief Energy
	mutable double energy;
};

#endif //LAGRANGIANHLLC_HPP
