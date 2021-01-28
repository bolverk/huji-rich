#include "hdf5_logger.hpp"
#include "../misc/utils.hpp"
#include "tessellation.hpp"
#include "../misc/hdf5_utils.hpp"
#include "../misc/serial_generate.hpp"
#include <numeric>

HDF5Logger::HDF5Logger(const string& fname):
  fname_(fname) {}

namespace{
  vector<int> sequential_integers(int n)
  {
    vector<int> res(static_cast<size_t>(n));
    iota(res.begin(), res.end(), 0);
    return res;
  }
}

void HDF5Logger::output(const Tessellation& v)
{
  const vector<Vector2D> points=
    serial_generate<int, Vector2D>
    (sequential_integers(v.GetPointNo()),
     [&v](int n){return v.GetMeshPoint(n);});
    
  HDF5Shortcut hdf5sc(fname_);
  hdf5sc("x", serial_generate<Vector2D, double>
	 (points,
	  [](const Vector2D& u){return u.x;}));
  hdf5sc("y", serial_generate<Vector2D, double>
	 (points,
	  [](const Vector2D& u){return u.y;}));
  //  hdf5sc("x",serial_generate(CoordinateExtractor(v,&Vector2D::x)));
  //  hdf5sc("y",serial_generate(CoordinateExtractor(v,&Vector2D::y)));
}
