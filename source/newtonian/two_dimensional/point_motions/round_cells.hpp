#ifndef ROUND_CELLS_HPP
#define ROUND_CELLS_HPP 1

#include "../point_motion.hpp"
#include "../../common/equation_of_state.hpp"

class RoundCells: public PointMotion
{
public:

  RoundCells(const PointMotion& pm,
	     const EquationOfState& eos,
	     double chi = 0.15,
	     double eta = 0.02);

  vector<Vector2D> operator()
  (const Tessellation& tess,
   const vector<ComputationalCell>& cells,
   double time) const;

private:

  Vector2D calc_dw(size_t i,
		   const Tessellation& tess,
		   const vector<ComputationalCell>& cells) const;

  const PointMotion& pm_;
  const EquationOfState& eos_;
  const double chi_;
  const double eta_;
};

#endif // ROUND_CELLS_HPP
