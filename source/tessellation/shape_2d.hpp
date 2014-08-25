/*! \file shape_2d.hpp
  \brief Two dimensional shapes
  \author Elad Steinberg
 */

#ifndef SHAPE_2D_HPP
#define SHAPE_2D_HPP 1

#include "geometry.hpp"

//! \brief Abstract type for a two dimensional shape
class Shape2D
{
public:

  /*! \brief Returns true is a point is inside the shape, false otherwise
    \param r Point in 2d space
    \return True if r is inside the shape
   */
  virtual bool operator()(Vector2D const& r) const = 0;

  //! \brief virtual destructor
  virtual ~Shape2D(void);
};

//! \brief A circle
class Circle: public Shape2D
{
public:

  /*! \brief Class constructor
    \param center Position of the center
    \param radius Radius
   */
  Circle(Vector2D const& center,
	 double radius);

  bool operator()(Vector2D const& r) const;

private:

  const Vector2D center_;
  const double radius_;
};

//! \brief Complement set of the points inside a certain shape
class Outside: public Shape2D
{
public:

  /*! \brief Class constructor
    \param shape Original shape
   */
  Outside(Shape2D const& shape);

  bool operator()(Vector2D const& r) const;

private:

  Shape2D const& shape_;
};

#endif // SHAPE_2D_HPP