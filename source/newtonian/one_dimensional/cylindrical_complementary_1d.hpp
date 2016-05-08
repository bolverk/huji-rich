#include "source_term_1d.hpp"

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
