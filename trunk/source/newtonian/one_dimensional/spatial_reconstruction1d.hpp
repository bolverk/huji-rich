/*! \file spatial_reconstruction1d.hpp
  \brief Spatial reconstruction
  \author Almog Yalinewich
 */

#ifndef SPATIAL_RECONSTRUCTION_1D_HPP
#define SPATIAL_RECONSTRUCTION_1D_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"

using std::vector;

//! \brief Base class for spatial reconstruction
class SpatialReconstruction1D
{
public:

  /*! \brief Returns the hydrodynamic state on the left or
  the right side of the interface position (=vertexes).
    \param vp Pointer to vertex positions
    \param hv Pointer to hydrodynamic variables
    \param i Vertex index
    \param interface_speed Speed of the interface between cells
    \param dir Direction (0 for left of the boundary, 1 for right ofthe boundary)
	\*** ( positive direction is from the left to the right).
    \param dt Time step
    \return Hydrodynamic variables on the left or right side of the interface
   */
  virtual Primitive InterpState(vector<double> const& vp, 
				vector<Primitive> const& hv, 
				double interface_speed,
				size_t i, int dir,double dt) const = 0;

  virtual ~SpatialReconstruction1D(void);
};

#endif // SPATIAL_RECONSTRUCTION_1D_HPP
