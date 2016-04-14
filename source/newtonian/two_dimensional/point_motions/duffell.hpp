#ifndef DUFFELL_HPP
#define DUFFELL_HPP 1

#include "../point_motion.hpp"

//! \brief Point motion based on Paul Duffell's scheme
//! \todo Add reference to paper
class Duffell: public PointMotion
{
public:

  /*! \brief Class constructor
    \param alpha Weight for average
    \param iter Number of iterations
   */
  Duffell
  (const double alpha,
   const int iter);

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time,TracerStickerNames const& tracerstickernames) const;

private:

  const double alpha_;
  const int iter_;
};

#endif // DUFFELL_HPP
