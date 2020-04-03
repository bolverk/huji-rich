#include "source_term_1d.hpp"

//! \brief Complementary source term for spherical simulations
//! \todo Change name to spherical complementary
class CylindricalComplementary1D: public SourceTerm1D
{
public:

  CylindricalComplementary1D(void);

  Conserved operator()
  (const SimulationState1D& state,
   size_t point,
   double t,
   double dt) const;
};
