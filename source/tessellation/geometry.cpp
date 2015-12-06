#include <cmath>
#include "geometry.hpp"

using namespace std;

Vector2D Rotate(Vector2D const& v, double a)
{
  return Vector2D(v.x*cos(a)-v.y*sin(a),
		  v.x*sin(a)+v.y*cos(a));
}

double Vector2D::distance(Vector2D const& v1) const
{
  return sqrt((x-v1.x)*(x-v1.x)
	      +(y-v1.y)*(y-v1.y));
}

double abs(Vector2D const& v)
{
  return sqrt(v.x*v.x+v.y*v.y);
}

Vector2D::Vector2D(void):
  x(0), y(0) {}

Vector2D::Vector2D(double ix, double iy):
  x(ix), y(iy) {}

Vector2D::Vector2D(const Vector2D& v):
  x(v.x), y(v.y) {}

void Vector2D::Set(double ix, double iy)
{
  x = ix;
  y = iy;
}

Vector2D& Vector2D::operator=(Vector2D const& v)
{
  if (this == &v)
    return *this;
  x=v.x;
  y=v.y;
  return *this;
}

Vector2D& Vector2D::operator*=(double s)
{
  x*=s;
  y*=s;
  return *this;
}

Vector2D& Vector2D::operator+=(Vector2D const& v)
{
  x += v.x;
  y += v.y;
  return *this;
}

Vector2D& Vector2D::operator-=(Vector2D const& v)
{
  x -= v.x;
  y -= v.y;
  return *this;
}

void Vector2D::Rotate(double a)
{
  x = x*cos(a)-y*sin(a);
  y = x*sin(a)+y*cos(a);
}

Vector2D operator+(Vector2D const& v1, Vector2D const& v2)
{
  return Vector2D(v1.x+v2.x,v1.y+v2.y);
}

Vector2D operator-(Vector2D const& v1,
		   Vector2D const& v2)
{
  return Vector2D(v1.x-v2.x,v1.y-v2.y);
}

Vector2D operator*(double d, Vector2D const& v)
{
  return Vector2D(v.x*d,v.y*d);
}

Vector2D operator*(Vector2D const& v, double d)
{
  return d*v;
}

Vector2D operator/(Vector2D const& v, double d)
{
  return Vector2D(v.x/d, v.y/d);
}

double ScalarProd(Vector2D const& v1,
		  Vector2D const& v2)
{
  return v1.x*v2.x +
    v1.y*v2.y;
}

double Projection(Vector2D const& v1, Vector2D const& v2)
{
  return ScalarProd(v1, v2)/abs(v2);
}

double CalcAngle(Vector2D const& v1, Vector2D const& v2)
{
  return acos(ScalarProd(v1, v2)/abs(v1)/abs(v2));
}

Vector2D Reflect(Vector2D const& v, Vector2D const& axis)
{
  return 2*ScalarProd(v, axis)*axis/pow(abs(axis),2)-v;
}

double distance(Vector2D const& v1, Vector2D const& v2)
{
  return abs(v1-v2);
}

double CrossProduct(Vector2D const& v1, Vector2D const& v2)
{
  return v1.x*v2.y-
    v1.y*v2.x;
}

Vector2D calc_mid_point(Vector2D const& v1, Vector2D const& v2)
{
  return Vector2D((v1.x+v2.x)*0.5, (v1.y+v2.y)*0.5);
}

Vector2D zcross(Vector2D const& v)
{
  return Vector2D(v.y,-v.x);
}

Vector2D pol2cart(double radius,
		  double angle)
{
  return Vector2D(radius*cos(angle),
		  radius*sin(angle));
}


Vector2D normalize(const Vector2D& v)
{
  return v/abs(v);
}

double dist_sqr(const Vector2D& v)
{
  return ScalarProd(v,v);
}

#ifdef RICH_MPI
size_t Vector2D::getChunkSize(void) const
{
  return 2;
}

vector<double> Vector2D::serialize(void) const
{
  vector<double> res(2,0);
  res.at(0) = x;
  res.at(1) = y;
  return res;
}

void Vector2D::unserialize
(const vector<double>& data)
{
  assert(data.size()==getChunkSize());
  x = data.at(0);
  y = data.at(1);
}
#endif // RICH_MPI
