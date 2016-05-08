#include "source_term_1d.hpp"

//! \brief Complementary source term for spherical simulations
//! \todo Change name to spherical complementary
class CylindricalComplementary1D: public ExternalForces1D
{
public:

  CylindricalComplementary1D(void);

  Conserved calc
  (const vector<double>& vertices,
   const vector<Primitive>& cells,
   int point,
   double t,
   double dt) const;
};
