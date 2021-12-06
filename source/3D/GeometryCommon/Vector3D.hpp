/*! \file Vector3D.hpp
\brief 3D Geometrical calculations
\author Itai Linial
*/

#ifndef Vector3D_HPP
#define Vector3D_HPP 1

#include <vector>
#include <limits>
#include <cmath>
#include "../../misc/serializable.hpp"

using std::vector;

//! \brief 3D Mathematical vector
class Vector3D : public Serializable
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
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D(double ix, double iy, double iz);

	/*! \brief Class copy constructor
	\param v Other vector
	*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D(const Vector3D& v);

	/*! \brief Set vector components
	\param ix x Component
	\param iy y Component
	\param iz z Component
	*/
	inline void Set(double ix, double iy, double iz) 
	{
		x = ix;
		y = iy;
		z = iz;
	}

	//! \brief Component in the x direction
	double x;

	//! \brief Component in the y direction
	double y;

	//! \brief Component in the z direction
	double z;

  /*! \brief Indexed access to member
    \param index Member index
    \return Reference to member
   */
	double& operator[](size_t index);

  /*! \brief Indexed access to member
    \param index Member index
    \return Value of member
   */
	double operator[](size_t index)const;

	/*! \brief Addition
	\param v Vector to be added
	\return Reference to sum
	*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D& operator+=(Vector3D const& v);

	/*! \brief Subtraction
	\param v Vector to be subtracted
	\return Difference
	*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D& operator-=(Vector3D const& v);

	/*! \brief Assignment operator
	\param v Vector to be copied
	\return The assigned value
	*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D& operator=(Vector3D const& v);
	
	/*! \brief Scalar product
	\param s Scalar
	\return Reference to the vector multiplied by scalar
	*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	Vector3D& operator*=(double s);

	/*! \brief Compare 3D-Vectors (up to an arbitrary precision)
	\param v Vector to be compared to
	\return True/False - according to the comparison results.
	*/
	bool operator==(Vector3D const& v) const;

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

  size_t getChunkSize(void) const override;
	
	vector<double> serialize(void) const override;

	void unserialize(const vector<double>& data) override;

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
	~Vector3D(void) override {}
};

/*! \brief Norm of a vector
\param v Three dimensional vector
\return Norm of v
*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
double abs(Vector3D const& v);

/*! \brief Norm of a vector, less accurate
\param v Three dimensional vector
\return Norm of v
*/
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
double fastabs(Vector3D const& v);

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
#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
inline double ScalarProd(Vector3D const& v1, Vector3D const& v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}


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
\param normal Normal to the reflection plane
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
inline Vector3D CrossProduct(Vector3D const& v1, Vector3D const& v2)
{
	return Vector3D(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

/*! \brief Cross product
  \param v1 First vector
  \param v2 Second vector
  \param res result
 */
inline void CrossProduct(Vector3D const& v1, Vector3D const& v2,Vector3D &res)
{
	res.x = v1.y*v2.z - v1.z*v2.y;
	res.y = v1.z*v2.x - v1.x*v2.z;
	res.z = v1.x*v2.y - v1.y*v2.x;
}


/*! \brief Splits a vector of 3D points to components
\param vIn Input vector of 3D points
\param vX Vector of x coordinates (out)
\param vY Vector of y coordinates (out)
\param vZ Vector of z coordinates (out)
*/
void Split(vector<Vector3D> const & vIn, vector<double> & vX, vector<double> & vY, vector<double> & vZ);

/*! \brief Normalise vector
  \param vec Vector
  \return Normalised vector
 */
Vector3D normalize(Vector3D const& vec);

#endif // Vector3D_HPP
