/*! \file PCM3D.hpp
\brief Class for pcm interpolation of the hydrodynamic variables
\author Elad Steinberg
*/

#ifndef PCM3D_HPP
#define PCM3D_HPP 1

#include "SpatialReconstruction3D.hpp"
#include "Ghost3D.hpp"

//! \brief Piecewise constant interpolation
class PCM3D : public SpatialReconstruction3D
{
private:
	Ghost3D const & ghost_;
	std::vector<Slope3D> slopes_;
public:
	/*! \brief Class constructor
	\param ghost The ghost point generator
	*/
	explicit PCM3D(Ghost3D const& ghost);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, double time,
		vector<pair<ComputationalCell3D, ComputationalCell3D> > &res, TracerStickerNames const& tracerstickersnames)const;

	void BuildSlopes(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, double time, TracerStickerNames const& tracerstickersnames);

	std::vector<Slope3D>& GetSlopes(void);
};

#endif // PCM3D_HPP
