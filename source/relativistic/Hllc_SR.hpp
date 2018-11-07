/*! \file Hllc_SR.hpp
\brief HLLC riemann solver for special relativity based on Mignone & Bodo 2005
\author Elad Steinberg
*/

#ifndef HLLCSR_HPP
#define HLLCSR_HPP 1

#include "../newtonian/common/riemann_solver.hpp"

//! \brief HLLC Riemann solver for an Eulerian grid
class Hllc_SR : public RiemannSolver
{
public:

	/*! \brief Class constructor
	 */
	explicit Hllc_SR();

	Conserved operator()
		(Primitive const& left,
			Primitive const& right,
			double velocity) const;
};

#endif //HLLCSR_HPP
