/*! \file pcm_scalar.hpp
  \brief Piecewise constant scalar interpolation
  \author Elad Steinberg
  \deprecated We are now using unified hydro + tracter interpolation
 */

#include "scalar_interpolation.hpp"

//! \brief PCM interpolation for tracers
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
