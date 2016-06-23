/*! \file point_motion_3d.hpp
  \brief Abstract class for the motion of the mesh generating points
  \author Almog Yalinewich
 */

#ifndef POINT_MOTION3D_HPP
#define POINT_MOTION3D_HPP 1

#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"

//! \brief Abstract class for point motion
class PointMotion3D
{
public:

	/*! \brief Calculates the velocity of all mesh points
	\param tess The tessellation
	\param cells Hydrodynamics cells
	\param time The simulation time
	\param tracerstickernames The names of the tracers and stickers
	\param res Velocities of the points, given as output
	*/
	virtual void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const = 0;

	/*! \brief Applies a small fix to the velocity of all mesh points once the time step is known
	\param tess The tessellation
	\param cells Hydrodynamics cells
	\param time The simulation time
	\param velocities Velocities of the points given as input and output
	\param dt The time step
	\param tracerstickernames The names of the tracers and stickers
	*/
	virtual void ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
		double dt, vector<Vector3D> &velocities, TracerStickerNames const& tracerstickernames)const;

  //! \brief Class destructor
  virtual ~PointMotion3D(void);
};

#endif // POINT_MOTION3D_HPP
