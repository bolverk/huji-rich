#include "hdf5_logger.hpp"
#include "../misc/utils.hpp"
#include "../misc/lazy_list.hpp"
#include "tessellation.hpp"
#include "../misc/hdf5_utils.hpp"

HDF5Logger::HDF5Logger(const string& fname):
  fname_(fname) {}

namespace{
  class CoordinateExtractor: public LazyList<double>
  {
  public:

    CoordinateExtractor(const Tessellation& tess,
			double Vector2D::* component):
      tess_(tess), component_(component) {}

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

    double operator[](size_t i) const
    {
      return (tess_.GetMeshPoint(static_cast<int>(i))).*component_;
    }

  private:
    const Tessellation& tess_;
    double Vector2D::* component_;
  };
}

void HDF5Logger::output(const VoronoiMesh& v)
{
  HDF5Shortcut hdf5sc(fname_);
  hdf5sc(string("x"),serial_generate(CoordinateExtractor(reinterpret_cast<const Tessellation&>(v),&Vector2D::x)));
  hdf5sc(string("y"),serial_generate(CoordinateExtractor(reinterpret_cast<const Tessellation&>(v),&Vector2D::y)));
}

void HDF5Logger::output(const Tessellation& v)
{
  HDF5Shortcut hdf5sc(fname_);
  hdf5sc("x",serial_generate(CoordinateExtractor(v,&Vector2D::x)));
  hdf5sc("y",serial_generate(CoordinateExtractor(v,&Vector2D::y)));
}
