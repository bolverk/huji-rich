/*! \file ersig.hpp
  \brief Exact Riemann solver for ideal gas equation of state
  \author Almog Yalinewich
*/

#ifndef ERSIG_HPP
#define ERSIG_HPP 1

#include <string>
#include "riemann_solver.hpp"

using std::string;

//! \brief Exact Riemann solver for ideal gas
class ERSIG: public RiemannSolver
{
public:

  /*! \brief Class constructor
    \param g Adiabatic index
    \param vacuum_behaviour Toggles optional behaviour in case of vacuum
   */
  ERSIG(double g, string const& vacuum_behaviour = "throw exception");

  Conserved operator()(Primitive const& left,
		       Primitive const& right,
		       double velocity) const;

private:

  const double g_;
  const string vacuum_behaviour_;
};

#endif // ERSIG_HPP
