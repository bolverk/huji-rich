/*! \file eulerian_3d.hpp
  \brief Eulerian point motion
  \author Almog Yalinewich
 */

#ifndef EULERIAN_3D_HPP
#define EULERIAN_3D_HPP 1

#include "point_motion_3d.hpp"

//! \brief Eulerian point motion
class Eulerian3D: public PointMotion3D
{
public:

  //! \brief Class constructor
  Eulerian3D(void);

  Vector3D operator()(const Vector3D& pos) const;
};

#endif // EULERIAN_3D_HPP
