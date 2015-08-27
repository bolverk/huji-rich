/*! \file LinearGaussImproved.hpp
\brief Linear interpolation that guarantees compliance with the equation of state and calcualtes the GG gradient from the CM of the cell
\author Elad Steinberg
*/

#ifndef LINEAR_GAUSS_IMPROVED
#define LINEAR_GAUSS_IMPROVED 1

#include "../../common/equation_of_state.hpp"
#include "../spatial_reconstruction.hpp"
#include <cmath>
#include "../../../misc/universal_error.hpp"

class LinearGaussImproved : public SpatialReconstruction
{
public:

	/*! \brief Class constructor
	\param eos Equation of state
	\param slf Slope limiter flag
	\param delta_v The GradV*L/Cs ratio needed for slope limiter
	\param theta The theta from tess in slope limiter.
	\param delta_P The pressure ratio for shock detection
	*/
	LinearGaussImproved(EquationOfState const& eos,bool slf = true,double delta_v = 0.2, 
		double theta = 0.5,double delta_P = 0.7);

	vector<pair<ComputationalCell, ComputationalCell> > operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells)const;

	/*! \brief Interpolates a cell
	\param cell The primitives of the cell
	\param cell_index The index of the cell
	\param cm The cell's center of mass
	\param target The location of the interpolation
	\return The interpolated value
	*/

	ComputationalCell Interp(ComputationalCell const& cell,size_t cell_index, Vector2D const& cm, Vector2D const& target)const;
private:
	EquationOfState const& eos_;
	mutable vector<pair<ComputationalCell, ComputationalCell> > rslopes_;
	bool slf_;
	double shockratio_, diffusecoeff_, pressure_ratio_;

	LinearGaussImproved(const LinearGaussImproved& origin);
	LinearGaussImproved& operator=(const LinearGaussImproved& origin);

};
#endif //LINEAR_GAUSS_IMPROVED
