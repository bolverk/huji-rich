/*! \file LagrangianHLLC.hpp
\brief HLLC riemann solver on an eulerian grid that assumes no mass transfer
\author Elad Steinberg
*/

#ifndef LAGRANGIANHLLC_HPP
#define LAGRANGIANHLLC_HPP 1

#include "riemann_solver.hpp"

//! \brief LagrangianHLLC Riemann solver for an Eulerian grid
class LagrangianHLLC : public RiemannSolver
{
public:

	Conserved operator()
		(Primitive const& left,
			Primitive const& right,
			double velocity) const;
};

#endif //LAGRANGIANHLLC_HPP
