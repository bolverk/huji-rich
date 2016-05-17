/*! \file zero_force_1d.hpp
  \brief Abstract class for zero external force
  \author Almog Yalinewich
*/

#ifndef ZERO_FORCE_1D_HPP
#define ZERO_FORCE_1D_HPP 1

#include "source_term_1d.hpp"

//! \brief Zero external force
class ZeroForce1D: public SourceTerm1D
{
public:
  Conserved operator()
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   size_t point,
   double t,
   double dt) const;
};

#endif // ZERO_FORCE_1D_HPP
