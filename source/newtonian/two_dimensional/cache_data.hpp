#ifndef CACHE_DATA_HPP
#define CACHE_DATA_HPP 1

#include "../../misc/cached_lazy_list.hpp"
#include "../../tessellation/tessellation.hpp"
#include "physical_geometry.hpp"

class CacheData
{
private:

  class VolumeCalculator: public LazyList<double>
  {
  public:

    VolumeCalculator(const Tessellation& tess,
		     const PhysicalGeometry& pg);

    size_t size(void) const;

    double operator[](const size_t i) const;

  private:
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };

  class AreaCalculator: public LazyList<double>
  {
  public:

    AreaCalculator(const Tessellation& tess,
		   const PhysicalGeometry& pg);

    size_t size(void) const;

    double operator[](const size_t i) const;

  private:
    const Tessellation& tess_;
    const PhysicalGeometry& pg_;
  };

  const VolumeCalculator volume_func_;
  const AreaCalculator area_func_;

public:

  CacheData(const Tessellation& tess,
	    const PhysicalGeometry& pg);

  void reset(void) const;

  const CachedLazyList<double> volumes;
  const CachedLazyList<double> areas;
};

#endif // CACHE_DATA_HPP
