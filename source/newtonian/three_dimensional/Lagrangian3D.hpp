/*! \file Lagrangian3D.hpp
\brief Lagrangian point motion scheme
\details Sets all the components velocities of all points to be the fluid's velocity
*/

#ifndef LAGRANGIAN3D_HPP
#define LAGRANGIAN3D_HPP 1

#include "point_motion_3d.hpp"

//! \brief Motion scheme where all point velocities are equal to the material velocity
class Lagrangian3D : public PointMotion3D
{
public:

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const;
};

#endif // LAGRANGIAN_HPP
