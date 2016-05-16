/*! \file round_cells.hpp
  \author Almog Yalinewich
  \brief Correction to point velocities that keeps cells round
 */

#ifndef ROUND_CELLS_HPP
#define ROUND_CELLS_HPP 1

#include "../point_motion.hpp"
#include "../../common/equation_of_state.hpp"
#include "../OuterBoundary.hpp"
#include "../geometric_outer_boundaries/PeriodicBox.hpp"

//! \brief Correction to point velocities that keeps cells round
//! \details Based on equation 63 in the Arepo paper
class RoundCells: public PointMotion
{
public:

  /*! \brief Class constructor
    \param pm Base point motion
    \param eos Equation of state
	\param outer The outer boudnary conditions, used for preventing points from getting outside the box. If periodic then no fix is applied
    \param chi chi parameter in equation 63
    \param eta eta parameter in equation 63
	\param cold Switch for cold flows
   */
	RoundCells(const PointMotion& pm, const EquationOfState& eos, OuterBoundary const& outer,
		double chi = 0.15, double eta = 0.02, bool cold = false);


	/*! \brief Class constructor
	\param pm Base point motion
	\param eos Equation of state
	\param chi chi parameter in equation 63
	\param eta eta parameter in equation 63
	\param cold Switch for cold flows
	*/
	RoundCells(const PointMotion& pm, const EquationOfState& eos, double chi = 0.15, double eta = 0.02,bool cold = false);

  vector<Vector2D> operator()(const Tessellation& tess,const vector<ComputationalCell>& cells,double time,
	  TracerStickerNames const& tracerstickernames) const;

  vector<Vector2D> ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	  double dt, vector<Vector2D> const& velocities, TracerStickerNames const& tracerstickernames)const;
private:

  Vector2D calc_dw(size_t i, const Tessellation& tess,const vector<ComputationalCell>& cells,
	  TracerStickerNames const& tracerstickernames) const;

  Vector2D calc_dw(size_t i, const Tessellation& tess, double dt,vector<ComputationalCell> const& cells,
	  TracerStickerNames const& tracerstickernames)const;

  const PointMotion& pm_;
  const EquationOfState& eos_;
  PeriodicBox pouter_;
  OuterBoundary const& outer_;
  const double chi_;
  const double eta_;
  const bool cold_;
};

#endif // ROUND_CELLS_HPP
