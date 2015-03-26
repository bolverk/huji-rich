#ifndef EXTENSIVE_UPDATER_HPP
#define EXTENSIVE_UPDATER_HPP 1

#include <vector>
#include "extensive.hpp"
#include "physical_geometry.hpp"
#include "../../tessellation/tessellation.hpp"
#include "cache_data.hpp"

using std::vector;

class ExtensiveUpdater
{
public:

  virtual void operator()
  (const vector<Extensive>& fluxes,
   const PhysicalGeometry& pg,
   const Tessellation& tess,
   const double dt,
   const CacheData& cd,
   vector<Extensive>& extensives) const = 0;

  virtual ~ExtensiveUpdater(void);
};
  
#endif // EXTENSIVE_UPDATER_HPP
