#ifndef DUFFELL_HPP
#define DUFFELL_HPP 1

#include "../point_motion.hpp"

class Duffell: public PointMotion
{
public:

  Duffell
  (const double alpha,
   const int iter);

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const;

private:

  const double alpha_;
  const int iter_;
};

#endif // DUFFELL_HPP
