/*! \file vertex_motion.hpp
  \brief Arbitrary motion of the vertices
  \author Almog Yalinewich
 */

#ifndef VERTEX_MOTION_HPP
#define VERTEX_MOTION_HPP 1

#include <vector>
#include "../common/hydrodynamic_variables.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"

using std::vector;

//! \brief Base class for vertex motion
class VertexMotion
{
public:
  
  /*! \brief Calculates the velocity of a vertex
    \param i Vertex index
    \param vp Pointer to vertices position container
    \param hv Pointer to hydrodynamics variables container
    \return Velocity of the vertex
  */
  virtual double operator()
  (int i, vector<double> const& vp,
   vector<ComputationalCell> const& hv) const = 0;

  virtual ~VertexMotion(void);
};

#endif
