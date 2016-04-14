/*! \file ConstantGravity.hpp
  \brief Constant gravity acceleration
  \author Elad Steinberg
*/

#ifndef CONSTANTGRAVITY_HPP
#define CONSTANTGRAVITY_HPP 1

#include "ConservativeForce.hpp"

//! \brief Acceleration due to constant gravity
class ConstantGravity: public Acceleration
{
public:
  /*! \brief Class constructor
    \param force The acceleration vector
  */
  explicit ConstantGravity(Vector2D const& force);
  
  Vector2D operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const double time,
   const int point,
   TracerStickerNames const& tracerstickernames) const;
 
private:
  const Vector2D force_;
};

#endif // CONSTANTGRAVITY_HPP
