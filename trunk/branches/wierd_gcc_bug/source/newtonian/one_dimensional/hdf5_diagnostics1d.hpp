#include <H5Cpp.h>
#include <string>
#include "hdsim.hpp"

using namespace std;
using namespace H5;

namespace diagnostics1d{
  void write_snapshot_to_hdf5(hdsim1D const& sim,
			      string const& fname);
}
