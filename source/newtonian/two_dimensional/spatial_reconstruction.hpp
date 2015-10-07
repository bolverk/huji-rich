/*! \file spatial_reconstruction.hpp
  \brief Abstract class for interpolation of the hydrodynamic variables
  \author Almog Yalinewich
*/

#ifndef SPATIAL_RECONSTRUCTION_HPP
#define SPATIAL_RECONSTRUCTION_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell_2d.hpp"

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
    \return List of pairs of primitive values on each cell
   */
  virtual vector<pair<ComputationalCell, ComputationalCell> >
  operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const = 0;

  virtual ~SpatialReconstruction(void);
};

#endif // SPATIAL_RECONSTRUCTION_HPP
