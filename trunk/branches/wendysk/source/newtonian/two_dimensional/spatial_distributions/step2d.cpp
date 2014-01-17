#include "step2d.hpp"

Step2D::Step2D(double xl, double xr,
	       double yd, double yu,
	       double vs, double vo):
  _xl(xl), _xr(xr), _yd(yd),
  _yu(yu), _vs(vs), _vo(vo) {}

double Step2D::EvalAt(Vector2D const& r) const
{
  if(r.x>_xl&&r.x<_xr&&
     r.y>_yd&&r.y<_yu)
    return _vs;
  else
    return _vo;
}
