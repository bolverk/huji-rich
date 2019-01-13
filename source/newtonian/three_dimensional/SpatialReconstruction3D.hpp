/*! \file SpatialReconstruction3D.hpp
\brief Abstract class for interpolation of the hydrodynamic variables
\author Elad Steinberg
*/

#ifndef SPATIAL_RECONSTRUCTION3D_HPP
#define SPATIAL_RECONSTRUCTION3D_HPP 1

#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "../three_dimensional/computational_cell.hpp"

using std::pair;

/*! \brief Spatial reconstruction of the primitive functions
\author Elad Steinberg
*/
class SpatialReconstruction3D
{
public:

	/*! \brief interpolates values on both sides of each interface
	\param tess Tessellation
	\param cells Computational cells
	\param time The sim time
	\param res List of pairs of primitive values on each edge given as output
	\param tracerstickersnames The names of the tracers and stickers
	*/
	virtual void operator()(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,double time,
		vector<pair<ComputationalCell3D, ComputationalCell3D> > &res,TracerStickerNames const& tracerstickersnames) const = 0;

	virtual ~SpatialReconstruction3D(void);

	virtual void BuildSlopes(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells,double time, TracerStickerNames const& tracerstickersnames) = 0;

	/*!
	\brief Returns the gradients
	\return The gradients
	*/
	virtual std::vector<Slope3D>& GetSlopes(void)=0;
};

#endif // SPATIAL_RECONSTRUCTION3D_HPP

