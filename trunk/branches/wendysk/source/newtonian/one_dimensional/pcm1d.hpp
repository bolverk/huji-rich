/*! \brief Piecewise constant interpolatio method for 1d simulations
  \author Almog Yalinewich
 */

#ifndef PCM1D_HPP
#define PCM1D_HPP 1

#include "spatial_reconstruction1d.hpp"

class PCM1D: public SpatialReconstruction1D
{
public:
  Primitive InterpState(vector<double> const& vp,
			vector<Primitive> const& hv,
			int i, int dir, double dt) const;
};

#endif // PCM1D_HPP 
