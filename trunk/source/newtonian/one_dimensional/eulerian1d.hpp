/*! \file eulerian1d.hpp
  \brief A velocity scheme where the velocities of all the vertices are 0
  \author Almog Yalinewich
*/

#ifndef EULERIAN1D_HPP
#define EULERIAN1D_HPP 1

#include "vertex_motion.hpp"

using std::vector;

//! \brief Eulerian vertex motion (no vertex motion)
class Eulerian1D: public VertexMotion
{
public:
  
  double CalcVelocity(int i, vector<double> const& vp,
		      vector<Primitive> const& hv) const;
};

#endif // EULERIAN1D_HPP
