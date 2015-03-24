/*! \file round_cells.hpp
  \author Almog Yalinewich
  \brief Correction to point velocities that keeps cells round
 */

#ifndef ROUND_CELLS_HPP
#define ROUND_CELLS_HPP 1

#include "../point_motion.hpp"
#include "../../common/equation_of_state.hpp"

//! \brief Correction to point velocities that keeps cells round
//! \details Based on equation 63 in the Arepo paper
class RoundCells: public PointMotion
{
public:

  /*! \brief Class constructor
    \param pm Base point motion
    \param eos Equation of state
    \param chi chi parameter in equation 63
    \param eta eta parameter in equation 63
   */
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
