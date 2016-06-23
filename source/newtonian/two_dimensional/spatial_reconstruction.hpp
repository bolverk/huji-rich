/*! \file spatial_reconstruction.hpp
  \brief Abstract class for interpolation of the hydrodynamic variables
  \author Almog Yalinewich
*/

#ifndef SPATIAL_RECONSTRUCTION_HPP
#define SPATIAL_RECONSTRUCTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell_2d.hpp"
#include "cache_data.hpp"

using std::pair;

/*! \brief Calculates the central point of the edge
  \param edge Edge
  \return Center of the edge
*/
Vector2D CalcCentroid(Edge const& edge);

/*! \brief Spatial reconstruction of the primitive functions
  \author Almog Yalinewich
*/
class SpatialReconstruction
{
public:

  /*! \brief interpolates values on both sides of each interface
    \param tess Tessellation
    \param cells Computational cells
	\param time The sim time
    \param res List of pairs of primitive values on each edge given as output
	\param tracerstickersnames The names of the tracers and stickers
	\param cd The cache data of the geometry
   */
  virtual void
  operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time, vector<pair<ComputationalCell, ComputationalCell> > &res,
	  TracerStickerNames const& tracerstickersnames,CacheData const& cd) const = 0;

  virtual ~SpatialReconstruction(void);
};

#endif // SPATIAL_RECONSTRUCTION_HPP
