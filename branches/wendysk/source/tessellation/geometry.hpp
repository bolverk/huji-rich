/*! \brief Geometrical calculations
  \author Almog Yalinewich
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1
#include <vector>
#include <boost/array.hpp>

//! \brief 2D Mathematical vector
class Vector2D
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

  /*! \brief Returns the x component 
    \returns The x component
  */
  double get_x(void) const;

  /*! \brief Returns the x component 
    \returns The x component
  */
  double get_y(void) const;

  /*! \brief Sets the X component 
    \param X The data to set
  */
  void set_x(double X);

  /*! \brief Sets the Y component 
  \param Y The data to set
  */
  void set_y(double Y);

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

  Vector2D& operator*=(double s);

  /*! \brief Rotates the vector in an anticlockwise direction
    \param a Angle of rotation (in radians)
   */
  void Rotate(double a);
  //! \brief Caluclates the distance from the Vector to v1 \param v1 The vector whose distance from is calculated \returns The distance
  double distance(Vector2D const& v1) const;
};

double abs(Vector2D const& v);

Vector2D operator+(Vector2D const& v1, Vector2D const& v2);

Vector2D operator-(Vector2D const& v1, Vector2D const& v2);

Vector2D operator*(double d, Vector2D const& v);

Vector2D operator*(Vector2D const& v, double d);

Vector2D operator/(Vector2D const& v, double d);

/*vector<double> operator*(vector<double> const& v,double c);

vector<double> operator*(double c,vector<double> const& v);

vector<double>& operator+=(vector<double> &lhs,vector<double> const&rhs);

vector<double>& operator-=(vector<double> &lhs,vector<double> const&rhs);
*/
double ScalarProd(Vector2D const& v1, Vector2D const& v2);

/*! \brief Returns the angle between two vectors (in radians)
 */
double CalcAngle(Vector2D const& v1, Vector2D const& v2);

double Projection(Vector2D const& v1, Vector2D const& v2);

Vector2D Rotate(Vector2D const& v, double a);

Vector2D Reflect(Vector2D const& v, Vector2D const& axis);

double distance(Vector2D const& v1, Vector2D const& v2);

/*! \brief Returns the z component of the cross product of two vectors
  \param v1 First vector
  \param v2 Second vector
 */
double CrossProduct(Vector2D const& v1, Vector2D const& v2);

Vector2D calc_mid_point(Vector2D const& v1, Vector2D const& v2);

Vector2D zcross(Vector2D const& v);

#endif // GEOMETRY_HPP
