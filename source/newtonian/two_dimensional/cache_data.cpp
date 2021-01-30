#include "cache_data.hpp"
#include "../../tessellation/ConvexHull.hpp"
#include "../../misc/lazy_list.hpp"
#include "../../misc/serial_generate.hpp"

CacheData::VolumeCalculator::VolumeCalculator
(const Tessellation& tess,
 const PhysicalGeometry& pg):
  tess_(tess), pg_(pg) {}

size_t CacheData::VolumeCalculator::size(void) const
{
  return static_cast<size_t>(tess_.GetPointNo());
}

double CacheData::VolumeCalculator::operator[](const size_t i) const
{
  return pg_.calcVolume
    (serial_generate<int, Edge>
     (tess_.GetCellEdges(static_cast<int>(i)),
      [&](int n)
      {return tess_.GetEdge(n);}));
      
    //    (serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))));
}

CacheData::AreaCalculator::AreaCalculator
(const Tessellation& tess,
 const PhysicalGeometry& pg):
  tess_(tess), pg_(pg) {}

size_t CacheData::AreaCalculator::size(void) const
{
  return tess_.getAllEdges().size();
}

double CacheData::AreaCalculator::operator[](const size_t i) const
{
  return pg_.calcArea(tess_.getAllEdges()[i]);
}

CacheData::CacheData(const Tessellation& tess,
		     const PhysicalGeometry& pg):
  volume_func_(tess,pg), area_func_(tess,pg),cm_func_(tess,pg),
  volumes(volume_func_), areas(area_func_),CMs(cm_func_) {}

void CacheData::reset(void) const
{
  volumes.reset();
  areas.reset();
  CMs.reset();
}

CacheData::CMCalculator::CMCalculator
(const Tessellation& tess,
	const PhysicalGeometry& pg) :
	tess_(tess), pg_(pg) {}

size_t CacheData::CMCalculator::size(void) const
{
	return static_cast<size_t>(tess_.GetPointNo());
}

Vector2D CacheData::CMCalculator::operator[](const size_t i) const
{
	vector<Vector2D> chull;
	ConvexHull(chull, tess_, static_cast<int>(i));
	return pg_.calcCentroid(chull);
}
