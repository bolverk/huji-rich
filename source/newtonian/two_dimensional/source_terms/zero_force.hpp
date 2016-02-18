/*! \file zero_force.hpp
  \brief Zero external force module
  \author Elad Steinberg
*/
#ifndef ZEROFORCE_HPP
#define ZEROFORCE_HPP 1

#include "../SourceTerm.hpp"
#include "ConservativeForce.hpp"

//! \brief Zero external force module
class ZeroForce: public SourceTerm
{
public:

  vector<Extensive> operator()
  (const Tessellation& tess,
   const PhysicalGeometry& pg,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const vector<Vector2D>& point_velocities,
   const double t,
	 TracerStickerNames const& tracerstickersnames) const;
};

#endif // ZEROFORCE_HPP
