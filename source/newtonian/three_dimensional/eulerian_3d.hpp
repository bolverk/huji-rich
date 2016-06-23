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

  void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	  double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const;
};

#endif // EULERIAN_3D_HPP
