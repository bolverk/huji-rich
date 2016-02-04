/*! \file simple_extensive_updater.hpp
  \author Almog Yalinewich
  \brief Simple extensive variable updater
 */

#ifndef SIMPLE_EXTENSIVE_UPDATER_HPP
#define SIMPLE_EXTENSIVE_UPDATER_HPP 1

#include "extensive_updater.hpp"

//! \brief Simple extensive variable updater
class SimpleExtensiveUpdater: public ExtensiveUpdater
{
public:
  
  void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& pg,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   const vector<ComputationalCell>& cells,
   vector<Extensive>& extensives,
	  double time) const;
};

#endif // SIMPLE_EXTENSIVE_UPDATER_HPP
