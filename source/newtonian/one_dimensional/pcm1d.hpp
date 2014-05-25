/*! \file pcm1d.hpp
  \brief Piecewise constant interpolation method for 1d simulations
  \author Almog Yalinewich
 */

#ifndef PCM1D_HPP
#define PCM1D_HPP 1

#include "spatial_reconstruction1d.hpp"

//! \brief Piecewise constant method for 1d spatial reconstruction
class PCM1D: public SpatialReconstruction1D
{
public:
  Primitive InterpState(vector<double> const& vp,
			vector<Primitive> const& hv,
			double interface_speed,
			size_t i, int dir, double dt) const;
};

#endif // PCM1D_HPP 
