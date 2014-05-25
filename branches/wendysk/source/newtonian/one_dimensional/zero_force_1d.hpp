/*! \brief Abstract class for zero external force
  \author Almog Yalinewich
*/

#ifndef ZERO_FORCE_1D_HPP
#define ZERO_FORCE_1D_HPP 1

#include "source_term_1d.hpp"

class ZeroForce1D: public ExternalForces1D
{
public:
  Conserved calc
  (vector<double> const& vertices,
   vector<Primitive> const& cells,
   int point,
   double t,
   double dt) const;
};

#endif // ZERO_FORCE_1D_HPP
