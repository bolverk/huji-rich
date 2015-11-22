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

//! \brief Correction to point velocities that keeps cells round
//! \details Based on Philip Mocz method
//! \todo Make work with MPI
class CentroidMotion : public PointMotion
{
public:

	/*! \brief Class constructor
	\param reduction_factor The factor to reduce the correction velocity (1/reduction_factor is the number of iterations to fix)
	\param outer The outer boudnary conditions, used for preventing points from getting outside the box. If periodic then no fix is applied
	\param niter The number of correction iterations to apply
	*/
	CentroidMotion(double reduction_factor,OuterBoundary const& outer, size_t niter = 2);

	vector<Vector2D> operator()(const Tessellation& tess, const vector<ComputationalCell>& cells, double time) const;

	vector<Vector2D> ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
		double dt, vector<Vector2D> const& velocities)const;
private:

	const double reduce_factor_;
	OuterBoundary const& outer_;
	const size_t niter_;
};

#endif // CENTROIDMOTION_HPP

