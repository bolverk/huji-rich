/*! \file CenterGravity.hpp
  \brief Point source gravity force
  \author Elad Steinberg
*/

#ifndef CENTERGRAVITY_HPP
#define CENTERGRAVITY_HPP 1

#include "ConservativeForce.hpp"

//! \brief Gravitational acceleration due to a pointlike mass
class CenterGravity: public Acceleration
{
public:
  /*! \brief Class constructor
    \param M The mass of the point source
    \param Rmin The softenning length
    \param center The location of the point source
  */
  CenterGravity
  (double M,double Rmin,const Vector2D& center);

  Vector2D operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   const vector<Extensive>& fluxes,
   const double time,
   const int point,
   TracerStickerNames const& tracerstickernames) const;

private:
  const double M_;
  const double Rmin_;
  const Vector2D center_;
};

#endif // CENTERGRAVITY_HPP
