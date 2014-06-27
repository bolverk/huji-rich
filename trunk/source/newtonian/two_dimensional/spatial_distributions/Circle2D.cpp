#include "Circle2D.hpp"


Circle2D::Circle2D(double xc,double yc,double R, double in,double out):
_R(R),_in(in),_out(out),center(xc,yc) {}


Circle2D::~Circle2D(void)
{
}

double Circle2D::operator()(Vector2D const& r) const
{
	if(r.distance(center)<_R)
		return _in;
	else
		return _out;
}
