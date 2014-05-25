#include "scalar_interpolation.hpp"

class PCMScalar: public ScalarInterpolation
{
public:

  void Prepare(Tessellation const* tess,
	       vector<vector<double> > const& tracers,
	       double dt,
	       double time);

  vector<double> Interpolate(Tessellation const* tess,
			       vector<vector<double> > const& tracers,
			       double dt,
			       Edge const& edge,
			       int side,
			       SInterpolationType interptype) const;
};
