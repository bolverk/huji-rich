/*! \file LagrangianHLLC.hpp
\brief HLLC riemann solver that keeps track of transfered internal energy
\author Elad Steinberg
*/

#ifndef LAGRANGIANHLLC_HPP
#define LAGRANGIANHLLC_HPP 1

#include "riemann_solver.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"

//! \brief LagrangianHLLC Riemann solver for an Eulerian grid
class LagrangianHLLC : public RiemannSolver
{
private:
	const bool massflux_;
public:
  /*! \brief Class constructor
    \param massflux Whether to apply correction for mass flux
   */
  explicit LagrangianHLLC(bool massflux=true);

	Conserved operator()
		(Primitive const& left,
			Primitive const& right,
			double velocity) const;

  //! \brief Energy
	mutable double energy;
};

#endif //LAGRANGIANHLLC_HPP
