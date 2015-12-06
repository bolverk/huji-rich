/*! \file geometry.hpp
  \brief Geometrical calculations
  \author Almog Yalinewich
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1
#include <vector>
#include <boost/array.hpp>
#ifdef RICH_MPI
#include "../misc/serializable.hpp"
#endif // RICH_MPI

//! \brief 2D Mathematical vector
class Vector2D
#ifdef RICH_MPI
  : public Serializable
#endif // RICH_MPI
{
public:

  /*! \brief Null constructor
    \details Sets all components to 0
   */
  Vector2D(void);

  /*! \brief Class constructor
    \param ix x Component
    \param iy y Component
   */
  Vector2D(double ix, double iy);

  /*! \brief Class copy constructor
    \param v Other vector
   */
  Vector2D(const Vector2D& v);

  /*! \brief Set vector components
    \param ix x Component
    \param iy y Component
   */
  void Set(double ix, double iy);

  //! \brief Component in the x direction
  double x;

  //! \brief Component in the y direction
  double y;

  /*! \brief Addition
    \param v Vector to be added
    \return Reference to sum
   */
  Vector2D& operator+=(Vector2D const& v);

  /*! \brief Subtraction
    \param v Vector to be subtracted
    \return Difference
   */
  Vector2D& operator-=(Vector2D const& v);

  /*! \brief Assigment operator
    \param v Vector to be copied
    \return The assigned value
   */
  Vector2D& operator=(Vector2D const& v);

  /*! \brief Scalar product
    \param s Scalar
    \return Reference to the vector multiplied by scalar
   */
  Vector2D& operator*=(double s);

  /*! \brief Rotates the vector in an anticlockwise direction
    \param a Angle of rotation (in radians)
   */
  void Rotate(double a);
  //! \brief Caluclates the distance from the Vector to v1 \param v1 The vector whose distance from is calculated \returns The distance
  double distance(Vector2D const& v1) const;

#ifdef RICH_MPI
  /*! \brief Serializer
    \param ar Archiver
    \param int Version
   */
  template<class Archive>
  void serialize
  (Archive& ar, 
   const unsigned int /*version*/)
  {
    ar & x;
    ar & y;
  }

  vector<double> serialize(void) const;

  size_t getChunkSize(void) const;

  void unserialize
  (const vector<double>& data);
#endif // RICH_MPI
};

/*! \brief Norm of a vector
  \param v Two dimensional vector
  \return Norm of v
 */
double abs(Vector2D const& v);

/*! \brief Term by term addition
  \param v1 First vector
  \param v2 Second vector
  \return Sum
 */
Vector2D operator+(Vector2D const& v1, Vector2D const& v2);

/*! \brief Term by term subtraction
  \param v1 First vector
  \param v2 Second vector
  \return Difference
 */
Vector2D operator-(Vector2D const& v1, Vector2D const& v2);

/*! \brief Scalar product
  \param d Scalar
  \param v Vector
  \return Two dimensional vector
 */
Vector2D operator*(double d, Vector2D const& v);

/*! \brief Scalar product
  \param v Vector
  \param d Scalar
  \return Two dimensional vector
 */
Vector2D operator*(Vector2D const& v, double d);

/*! \brief Scalar division
  \param v Vector
  \param d Scalar
  \return Two dimensional vector
 */
Vector2D operator/(Vector2D const& v, double d);

/*! \brief Scalar product of two vectors
  \param v1 2D vector
  \param v2 2D vector
  \return Scalar product of v1 and v2
 */
double ScalarProd(Vector2D const& v1, Vector2D const& v2);

/*! \brief Returns the angle between two vectors (in radians)
  \param v1 First vector
  \param v2 Second vector
  \return Angle
 */
double CalcAngle(Vector2D const& v1, Vector2D const& v2);

/*! \brief Calculates the projection of one vector in the direction of the second
  \param v1 First vector
  \param v2 Direction of the projection
  \return Component of v1 in the direction of v2
 */
double Projection(Vector2D const& v1, Vector2D const& v2);

/*! \brief Rotates a vector
  \param v Vector
  \param a Angle
  \return Rotated vector
 */
Vector2D Rotate(Vector2D const& v, double a);

/*! \brief Reflect vector
  \param v Vector
  \param axis Axis of reflection
  \return Reflection of v about axis
 */
Vector2D Reflect(Vector2D const& v, Vector2D const& axis);

/*! \brief Calculates the distance between two vectors
  \param v1 First vector
  \param v2 Second vector
  \return distance between v1 and v2
 */
double distance(Vector2D const& v1, Vector2D const& v2);

/*! \brief Returns the z component of the cross product of two vectors
  \param v1 First vector
  \param v2 Second vector
  \return z component of the cross product between v1 and v2
 */
double CrossProduct(Vector2D const& v1, Vector2D const& v2);

/*! \brief Calculates the mid point between two vectors
  \param v1 First vector
  \param v2 Second vector
  \return Distance between v1 and v2
 */
Vector2D calc_mid_point(Vector2D const& v1, Vector2D const& v2);

/*! \brief Cross product of a vector in x,y plane with a unit vector in the z direction
  \param v Vector in the x,y plane
  \return Two dimensional vector
 */
Vector2D zcross(Vector2D const& v);

/*! \brief Converts from polar coordinates to cartesian coordinates
  \param radius Radius
  \param angle Angle relative to the x axis
  \return Same vector in cartesian coordiantes
 */
Vector2D pol2cart(double radius, double angle);

/*! \brief Normalized a vector
  \param v Original vector
  \return Vector divided by its norm
 */
Vector2D normalize(const Vector2D& v);

/*! \brief Calculates the square of the distance. This is computationaly cheaper then actually calculating the distance
  \param v Vector
  \return Square of the distance
 */
double dist_sqr(const Vector2D& v);

#endif // GEOMETRY_HPP
