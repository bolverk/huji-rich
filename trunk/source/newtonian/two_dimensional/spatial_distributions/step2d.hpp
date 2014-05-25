/*! \file step2d.hpp
  \brief A spatial distribution with a square step function
  \author Almog Yalinewich
 */

#ifndef STEP2D_HPP
#define STEP2D_HPP 1

#include "../spatial_distribution2d.hpp"

//! \brief A spatial distribution with a square step function
class Step2D: public SpatialDistribution
{
private:

  double _xl, _xr, _yd, _yu, _vs, _vo;

public:

  /*! \brief Class constructor
    \param xl Left boundary of the step
    \param xr Right boundary of the step
    \param yd Lower boundary of the step
    \param yu Upper boundary of the step
    \param vs Value on the step
    \param vo Value outside the step
   */
  Step2D(double xl, double xr,
	 double yd, double yu,
	 double vs, double vo);

  double EvalAt(Vector2D const& r) const;
};

#endif // STEP2D_HPP
