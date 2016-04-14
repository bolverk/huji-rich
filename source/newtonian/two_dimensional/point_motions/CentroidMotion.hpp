/*! \file CentroidMotion.hpp
\author Elad Steinberg
\brief Correction to point velocities that keeps cells round by prediction of cell center
*/

#ifndef CENTROIDMOTION_HPP
#define CENTROIDMOTION_HPP 1

#include "../point_motion.hpp"
#include "../../common/equation_of_state.hpp"
#include "../OuterBoundary.hpp"
#include "../geometric_outer_boundaries/PeriodicBox.hpp"
#include "../interpolations/LinearGaussImproved.hpp"

//! \brief Correction to point velocities that keeps cells round
//! \details Based on Philip Mocz method
//! \todo Make work with MPI
class CentroidMotion : public PointMotion
{
public:

	/*! \brief Class constructor
	\param reduction_factor The factor to reduce the correction velocity (1/reduction_factor is the number of iterations to fix)
	\param niter The number of correction iterations to apply
	\param coldflow Flag if to use the soundspeed for the correction or the time step if the flow is cold
	\param eos The equation of state
	\param toignore List of sticker names not to apply correction for
	*/
	CentroidMotion(double reduction_factor,EquationOfState const& eos, bool coldflow = false, size_t niter = 2,
		vector<string> toignore = vector<string>());

	vector<Vector2D> operator()(const Tessellation& tess, const vector<ComputationalCell>& cells, double time,
		TracerStickerNames const& tracerstickernames) const;

	vector<Vector2D> ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
		double dt, vector<Vector2D> const& velocities, TracerStickerNames const& tracerstickernames)const;
private:

	const double reduce_factor_;
	EquationOfState const& eos_;
	const bool cold_;
	const size_t niter_;
	const vector<string> toignore_;
};

#endif // CENTROIDMOTION_HPP

