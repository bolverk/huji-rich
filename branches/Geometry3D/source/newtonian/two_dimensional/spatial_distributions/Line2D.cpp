#include "Line2D.hpp"

Line2D::Line2D(double a,double b,double above,double under):
_a(a),_b(b),_above(above),_under(under){}

Line2D::~Line2D(void)
{
}

double Line2D::operator()(Vector2D const& r) const
{
	if(r.y>(_a*r.x+_b))
		return _above;
	else
		return _under;
}
