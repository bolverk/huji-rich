/*! \file pcm2d.hpp
  \brief A 2D piecewise contant interpolation method
  \author Almog Yalinewich
 */

#ifndef PCM2D_HPP
#define PCM2D_HPP 1

#include "../spatial_reconstruction.hpp"

//! \brief Piecewise constant method
class PCM2D: public SpatialReconstruction
{
public:

  void Prepare(Tessellation const& tessellation,
	       vector<Primitive> const& cells,
	       vector<vector<double> > const& tracers,
	       double dt,double time);

  Primitive Interpolate(Tessellation const& /*tessellation*/,
			vector<Primitive> const& cells,
			double /*dt*/,Edge const& edge, int side,
			InterpolationType /*interptype*/,Vector2D const& /*vface*/) const;

  vector<double> interpolateTracers
  (Tessellation const& tess,vector<Primitive> const& cells,
   vector<vector<double> > const& tracers,
   double dt,
   Edge const& edge,
   int side,
   InterpolationType interp_type,Vector2D const& vface) const;
};

#endif // PCM2D_HPP
