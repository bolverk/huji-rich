/*! \file hllc.hpp
  \brief HLLC riemann solver on an eulerian grid
  \details This file is based on a code originally written by Omer Bromberg
  \author Almog Yalinewich
*/

#ifndef HLLC_HPP
#define HLLC_HPP 1

#include "riemann_solver.hpp"

//! \brief HLLC Riemann solver for an Eulerian grid
class Hllc: public RiemannSolver
{
public:

  Conserved operator()
  (Primitive const& left,
   Primitive const& right,
   double velocity) const;
};

#endif //HLLC_EULERIAN_HPP
