/*! \brief A velocity scheme where the velocities of all vertices are determined by the velocities of adjacents cells
  \author Almog Yalinewich
*/

#ifndef LAGRANGIAN1D_HPP
#define LAGRANGIAN1D_HPP 1

#include "vertex_motion.hpp"

class Lagrangian1D: public VertexMotion
{
public:

  /*! \brief Class constructor
    \param rigid_walls Toggles rigid walls at the edges
   */
  Lagrangian1D(bool rigid_walls);

  double CalcVelocity(int i, vector<double> const& vp,
		      vector<Primitive> const& hv) const;

private:

  bool rigid_walls_;
};

#endif // LAGRANGIAN1D_HPP
