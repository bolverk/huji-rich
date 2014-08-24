#include "triangle_area.hpp"

double calc_triangle_area(const Vector2D& p1,
			  const Vector2D& p2,
			  const Vector2D& p3)
{
  return 0.5*abs(ScalarProd(p2-p1,zcross(p3-p1)));
}
