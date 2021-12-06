#include <cmath>
#include <math.h>
#include "Vector3D.hpp"
#include "../../misc/utils.hpp"

#define EPSILON 1e-12

namespace
{
	static inline double my_round(double val)
	{    
		return floor(val + 0.5);
	}
}

Vector3D RotateX(Vector3D const& v, double a)
{
	Vector3D res;
	res.x = v.x;
	res.y = v.y*cos(a) - v.z*sin(a);
	res.z = v.y*sin(a) + v.z*cos(a);
	return res;
}

Vector3D RotateY(Vector3D const& v, double a)
{
	Vector3D res;
	res.x = v.x*cos(a) + v.z*sin(a);
	res.y = v.y;
	res.z = -v.x*sin(a) + v.z*cos(a);
	return res;
}

Vector3D RotateZ(Vector3D const& v, double a)
{
	Vector3D res;
	res.x = v.x*cos(a) - v.y*sin(a);
	res.y = v.x*sin(a) + v.y*cos(a);
	res.z = v.z;
	return res;
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
double abs(Vector3D const& v)
{
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
double fastabs(Vector3D const& v)
{
	return fastsqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}


Vector3D::Vector3D(void) :
x(0), y(0), z(0) {}

Vector3D::Vector3D(double ix, double iy, double iz) :
x(ix), y(iy), z(iz) {}

Vector3D::Vector3D(const Vector3D& v) :
x(v.x), y(v.y), z(v.z) {}

size_t Vector3D::getChunkSize(void) const
{
	return 3;
}

vector<double> Vector3D::serialize(void) const
{
	vector<double> res(3);
	res[0] = x;
	res[1] = y;
	res[2] = z;
	return res;
}

void Vector3D::unserialize(const vector<double>& data)
{
	assert(data.size() == 3);
	x = data[0];
	y = data[1];
	z = data[2];
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
Vector3D& Vector3D::operator=(Vector3D const& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
	return *this;
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
Vector3D& Vector3D::operator*=(double s)
{
	x *= s;
	y *= s;
	z *= s;
	return *this;
}

double& Vector3D::operator[](size_t index)
{
	if (index == 0)
		return x;
	if (index == 1)
		return y;
	if (index == 2)
		return z;
	assert(false);
	double *temp = new double(std::numeric_limits<double>::signaling_NaN());
	return *temp;
}

double Vector3D::operator[](size_t index)const
{
	if (index == 0)
		return x;
	if (index == 1)
		return y;
	if (index == 2)
		return z;
	assert(false);
	return std::numeric_limits<double>::signaling_NaN();
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
Vector3D& Vector3D::operator+=(Vector3D const& v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
Vector3D& Vector3D::operator-=(Vector3D const& v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

// Note - since working with double precision, two vectors are assumed to be "equal",
// if their coordinates agree up to precision EPSILON
bool Vector3D::operator==(Vector3D const& v) const
{
	return (std::abs(x - v.x) < EPSILON) && (std::abs(y - v.y) < EPSILON) && (std::abs(z - v.z) < EPSILON);
}

void Vector3D::RotateX(double a)
{
	Vector3D v;
	v.x = x;
	v.y = y*cos(a) - z*sin(a);
	v.z = y*sin(a) + z*cos(a);

	*this = v;
}

void Vector3D::RotateY(double a)
{
	Vector3D v;
	v.x = x*cos(a) + z*sin(a);
	v.y = y;
	v.z = -x*sin(a) + z*cos(a);

	*this = v;
}

void Vector3D::RotateZ(double a)
{
	Vector3D v;
	v.x = x*cos(a) - y*sin(a);
	v.y = x*sin(a) + y*cos(a);
	v.z = z;

	*this = v;
}

void Vector3D::Round()
{
	x = my_round(x);
	y = my_round(y);
	z = my_round(z);
}

Vector3D operator+(Vector3D const& v1, Vector3D const& v2)
{
	Vector3D res;
	res.x = v1.x + v2.x;
	res.y = v1.y + v2.y;
	res.z = v1.z + v2.z;
	return res;
}

Vector3D operator-(Vector3D const& v1,
	Vector3D const& v2)
{
	Vector3D res;
	res.x = v1.x - v2.x;
	res.y = v1.y - v2.y;
	res.z = v1.z - v2.z;
	return res;
}

Vector3D operator*(double d, Vector3D const& v)
{
	Vector3D res;
	res.x = v.x * d;
	res.y = v.y * d;
	res.z = v.z * d;
	return res;
}

Vector3D operator*(Vector3D const& v, double d)
{
	return d*v;
}

Vector3D operator/(Vector3D const& v, double d)
{
	Vector3D res;
	res.x = v.x / d;
	res.y = v.y / d;
	res.z = v.z / d;
	return res;
}

double Projection(Vector3D const& v1, Vector3D const& v2)
{
	return ScalarProd(v1, v2) / abs(v2);
}

double CalcAngle(Vector3D const& v1, Vector3D const& v2)
{
	return acos(ScalarProd(v1, v2) / abs(v1) / abs(v2));
}

Vector3D Reflect(Vector3D const& v, Vector3D const& normal)
{
	return v - 2 * ScalarProd(v, normal)*normal / ScalarProd(normal,normal);
}

double distance(Vector3D const& v1, Vector3D const& v2)
{
	return abs(v1 - v2);
}

void Split(vector<Vector3D> const & vIn, vector<double> & vX, vector<double> & vY, vector<double> & vZ)
{
	vX.resize(vIn.size());
	vY.resize(vIn.size());
	vZ.resize(vIn.size());

	for (std::size_t ii = 0; ii < vIn.size(); ++ii)
	{
		vX[ii] = vIn[ii].x;
		vY[ii] = vIn[ii].y;
		vZ[ii] = vIn[ii].z;
	}
	return;
}

Vector3D normalize(Vector3D const& vec)
{
	double l = abs(vec);
	return vec / l;
}
