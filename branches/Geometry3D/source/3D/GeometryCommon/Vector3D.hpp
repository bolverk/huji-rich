/*! \file Vector3D.hpp
\brief 3D Geometrical calculations
\author Itai Linial
*/

#ifndef Vector3D_HPP
#define Vector3D_HPP 1

#include <vector>
#include <cmath>
using namespace std;

//! \brief 3D Mathematical vector
class Vector3D
{
public:

	/*! \brief Null constructor
	\details Sets all components to 0
	*/
	Vector3D(void);

	/*! \brief Class constructor
	\param ix x Component
	\param iy y Component
	\param iz z Component
	*/
	Vector3D(double ix, double iy, double iz);

	/*! \brief Class copy constructor
	\param v Other vector
	*/
	Vector3D(const Vector3D& v);

	/*! \brief Set vector components
	\param ix x Component
	\param iy y Component
	\param iz z Component
	*/
	void Set(double ix, double iy, double iz);

	//! \brief Component in the x direction
	double x;

	//! \brief Component in the y direction
	double y;

	//! \brief Component in the z direction
	double z;

	/*! \brief Addition
	\param v Vector to be added
	\return Reference to sum
	*/
	Vector3D& operator+=(Vector3D const& v);

	/*! \brief Subtraction
	\param v Vector to be subtracted
	\return Difference
	*/
	Vector3D& operator-=(Vector3D const& v);

	/*! \brief Assignment operator
	\param v Vector to be copied
	\return The assigned value
	*/
	Vector3D& operator=(Vector3D const& v);

	/*! \brief Scalar product
	\param s Scalar
	\return Reference to the vector multiplied by scalar
	*/
	Vector3D& operator*=(double s);

	/*! \brief Compare 3D-Vectors (up to an arbitrary precision)
	\param v Vector to be compared to
	\return True/False - according to the comparison results.
	*/
	bool operator==(const Vector3D &v) const;

	/*! \brief Rotates the vector around the X axes
	\param a Angle of rotation (in radians)
	*/
	void RotateX(double a);

	/*! \brief Rotates the vector around the Y axes
	\param a Angle of rotation (in radians)
	*/
	void RotateY(double a);

	/*! \brief Rotates the vector around the Z axes
	\param a Angle of rotation (in radians)
	*/
	void RotateZ(double a);

	/*! \brief Integer round of the vector's entries
	*/
	void Round();
};

/*! \brief Norm of a vector
\param v Three dimensional vector
\return Norm of v
*/
double abs(Vector3D const& v);

/*! \brief Term by term addition
\param v1 First vector
\param v2 Second vector
\return Sum
*/
Vector3D operator+(Vector3D const& v1, Vector3D const& v2);

/*! \brief Term by term subtraction
\param v1 First vector
\param v2 Second vector
\return Difference
*/
Vector3D operator-(Vector3D const& v1, Vector3D const& v2);

/*! \brief Scalar product
\param d Scalar
\param v Vector
\return Three dimensional vector
*/
Vector3D operator*(double d, Vector3D const& v);

/*! \brief Scalar product
\param v Vector
\param d Scalar
\return Three dimensional vector
*/
Vector3D operator*(Vector3D const& v, double d);

/*! \brief Scalar division
\param v Vector
\param d Scalar
\return Three dimensional vector
*/
Vector3D operator/(Vector3D const& v, double d);

/*! \brief Scalar product of two vectors
\param v1 3D vector
\param v2 3D vector
\return Scalar product of v1 and v2
*/
double ScalarProd(Vector3D const& v1, Vector3D const& v2);

/*! \brief Returns the angle between two vectors (in radians)
\param v1 First vector
\param v2 Second vector
\return Angle (radians)
*/
double CalcAngle(Vector3D const& v1, Vector3D const& v2);

/*! \brief Calculates the projection of one vector in the direction of the second
\param v1 First vector
\param v2 Direction of the projection
\return Component of v1 in the direction of v2
*/
double Projection(Vector3D const& v1, Vector3D const& v2);

/*! \brief Rotates a 3D-vector around the X axis
\param v Vector
\param a  (in radians)
\return Rotated vector
*/
Vector3D RotateX(Vector3D const& v, double a );

/*! \brief Rotates a 3D-vector around the Y axis
\param v Vector
\param a  (in radians)
\return Rotated vector
*/
Vector3D RotateY(Vector3D const& v, double a);

/*! \brief Rotates a 3D-vector around the Z axis
\param v Vector
\param a  (in radians)
\return Rotated vector
*/
Vector3D RotateZ(Vector3D const& v, double a);

/*! \brief Reflect vector
\param v Vector
\param axis Normal to the reflection plane
\return Reflection of v about axis
*/
Vector3D Reflect(Vector3D const& v, Vector3D const& normal);

/*! \brief Calculates the distance between two vectors
\param v1 First vector
\param v2 Second vector
\return distance between v1 and v2
*/
double distance(Vector3D const& v1, Vector3D const& v2);

/*! \brief Returns the cross product of two vectors
\param v1 First vector
\param v2 Second vector
\return Cross product between v1 and v2
*/
Vector3D CrossProduct(Vector3D const& v1, Vector3D const& v2);

/*! \brief Splits a vector of 3D points to components
\param vIn Input vector of 3D points
\param vX Vector of x coordinates (out)
\param vY Vector of y coordinates (out)
\param vZ Vector of z coordinates (out)
*/
void Split(vector<Vector3D> const & vIn, vector<double> & vX, vector<double> & vY, vector<double> & vZ);
#endif // Vector3D_HPP