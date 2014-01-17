/*! \brief Unifrom spatial distribution
  \author Almog Yalinewich
 */

#ifndef UNIFORM2D_HPP
#define UNIFORM2D_HPP 1

#include "../spatial_distribution2d.hpp"

//! \brief Uniform distribution
class Uniform2D: public SpatialDistribution
{
private:

  double _val;

public:

  /*! \brief Class constructr
    \param val Value
   */
  Uniform2D(double val);

  double EvalAt(Vector2D const& /*r*/) const;
};

#endif // UNIFORM2D_HPP
