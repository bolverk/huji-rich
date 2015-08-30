/*! \file arepo_interp.hpp
  \brief Interpolation described  in arepo paper
  \author Almog Yalinewich
 */

#ifndef AREPO_INTERP_HPP
#define AREPO_INTERP_HPP 1

#include "../common/ideal_gas.hpp"
#include "spatial_reconstruction1d.hpp"

//! \brief Interpolation based on Arepo's method
class ArepoInterp: public SpatialReconstruction1D
{
public:

  /*! \brief Class constructor
    \param eos Equation of state
   */
  explicit ArepoInterp(IdealGas const& eos);

  Primitive InterpState
  (vector<double> const& vp,
   vector<Primitive> const& hv,
   double interface_speed,
   size_t i, int dir, double dt) const;

private:
  IdealGas const& eos_;
};

#endif // AREPO_INTERP_HPP
