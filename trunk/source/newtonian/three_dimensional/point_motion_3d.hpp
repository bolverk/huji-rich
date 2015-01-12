/*! \file point_motion_3d.hpp
  \brief Abstract class for the motion of the mesh generating points
  \author Almog Yalinewich
 */

#ifndef POINT_MOTION3D_HPP
#define POINT_MOTION3D_HPP 1

#include "../../3D/GeometryCommon/Vector3D.hpp"

//! \brief Abstract class for point motion
class PointMotion3D
{
public:

  /*! \brief Returns the velocity of a mesh generating point
    \param pos Current position of mesh generating point
    \return Velocity of said point
   */
  virtual Vector3D operator()(const Vector3D& pos) const = 0;

  //! \brief Class destructor
  virtual ~PointMotion3D(void);
};

#endif // POINT_MOTION3D_HPP
