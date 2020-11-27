#include "source_term_1d.hpp"

//! \brief Complementary source term for spherical simulations
//! \todo Change name to spherical complementary
class CylindricalComplementary1D: public SourceTerm1D
{
public:

  CylindricalComplementary1D(void);

  Extensive operator()
  (const SimulationState1D& state,
   size_t point,
   const vector<Extensive>& fluxes,
   const PhysicalGeometry1D& pg,
   double t,
   double dt) const;
};
