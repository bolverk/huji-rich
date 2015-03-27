#include "cache_data.hpp"
#include "../../misc/lazy_list.hpp"

namespace {
  class CellEdgesGetter: public LazyList<Edge>
  {
  public:

    CellEdgesGetter(const Tessellation& tess, int n):
      tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

    size_t size(void) const
    {
      return edge_indices_.size();
    }

    Edge operator[](size_t i) const
    {
      return tess_.GetEdge(edge_indices_[i]);
    }

  private:
    const Tessellation& tess_;
    const vector<int> edge_indices_;
  };
}

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
    (serial_generate(CellEdgesGetter(tess_,static_cast<int>(i))));
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
  volume_func_(tess,pg), area_func_(tess,pg),
  volumes(volume_func_), areas(area_func_) {}

void CacheData::reset(void) const
{
  volumes.reset();
  areas.reset();
}
