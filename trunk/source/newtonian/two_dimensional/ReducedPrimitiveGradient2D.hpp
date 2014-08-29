/*! \file ReducedPrimitiveGradient2D.hpp
  \brief A 2D gradient of the hydrodynamic variables
  \author Almog Yalinewich
 */

#ifndef REDUCEDPRIMITIVEGRADIENT_HPP
#define REDUCEDPRIMITIVEGRADIENT_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include <vector>
#include "../../tessellation/geometry.hpp"

using namespace std;

//! \brief Gradient of the hydrodynamic variables, without energy and sound speed
class ReducedPrimitiveGradient2D
{
public:

  //! \brief Null constructor (sets everything to zero)
  ReducedPrimitiveGradient2D(void);

  /*! \brief Class constructor
    \param d Density gradient
    \param p Pressure gradient
    \param vx x velocity gradient
    \param vy y Velocity gradient
    \param trace Tracers gradient
   */
  ReducedPrimitiveGradient2D
  (Vector2D const& d,
   Vector2D const& p,
   Vector2D const& vx,
   Vector2D const& vy,
   vector<Vector2D> const& trace);

  //! \brief Density gradient
  Vector2D density;

  //! \brief Pressure gradient
  Vector2D pressure;

  //! \brief Gradient of the x velocity
  Vector2D xvelocity;

  //! \brief Gradient of the y velocity
  Vector2D yvelocity;

  //! \brief Gradient of the tracers
  vector<Vector2D> tracers;

  /*! \brief Term by term sum
    \param source Gradient to add
    \return sum
   */
  ReducedPrimitiveGradient2D& operator+=
  (ReducedPrimitiveGradient2D const& source);

  /*! \brief Scalar multiplication
    \param s Scalar
    \return Gradient multiplied by a scalar
   */
  ReducedPrimitiveGradient2D& operator*=
  (double s);
};

#endif //REDUCEDPRIMITIVEGRADIENT_HPP
