/*! \file LinearGauss3D.hpp
\brief Linear interpolation that guarantees compliance with the equation of state and calcualtes the GG gradient from the CM of the cell
\author Elad Steinberg
*/

#ifndef LINEAR_GAUSS3D_HPP
#define LINEAR_GAUSS3D_HPP 1

#include "../common/equation_of_state.hpp"
#include "SpatialReconstruction3D.hpp"
#include <cmath>
#include "../../misc/universal_error.hpp"
#include "Ghost3D.hpp"

//! \brief Linear gauss interpolation
class LinearGauss3D : public SpatialReconstruction3D
{
public:

	/*! \brief Class constructor
	\param eos Equation of state
	\param slf Slope limiter flag
	\param delta_v The GradV*L/Cs ratio needed for slope limiter
	\param theta The theta from tess in slope limiter.
	\param delta_P The pressure ratio for shock detection
	\param ghost The ghost point generator
	\param SR Flag for SR
	\param calc_tracers Names of tracers for which to calc the slope for
	\param skip_key The sticker name to skip cells for taking them into account for the slope limit
	\param tsn The names of the stickers and tracers
	*/
	LinearGauss3D(EquationOfState const& eos,TracerStickerNames const& tsn, Ghost3D const& ghost,bool slf = true,double delta_v = 0.2,
		double theta = 0.5,double delta_P = 0.7,bool SR=false,const vector<string>& calc_tracers = vector<string>(),
		string skip_key = string(),bool pressure_calc=true);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, double time,
		vector<pair<ComputationalCell3D, ComputationalCell3D> > &res, TracerStickerNames const& tracerstickersnames) const;

	/*! \brief Interpolates a cell
	\param cell The primitives of the cell
	\param cell_index The index of the cell
	\param cm The cell's center of mass
	\param target The location of the interpolation
	\param res The interpolated value
	\param tsn The names of the stickers and tracers
	\param eos The equation of state
	*/
	void Interp(ComputationalCell3D &res,ComputationalCell3D const& cell, size_t cell_index, Vector3D const& cm, Vector3D const& target,
		EquationOfState const& eos,TracerStickerNames const& tsn)const;

	/*!
	\brief Returns the gradients
	\return The gradients
	*/
	vector<Slope3D>& GetSlopes(void);

	/*!
	\brief Returns the unsloped limtied gradients
	\return The gradients
	*/
	vector<Slope3D>& GetSlopesUnlimited(void)const;

	void BuildSlopes(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, double time, TracerStickerNames const& tracerstickersnames);

private:
	EquationOfState const& eos_;
	TracerStickerNames const& tsn_;
	Ghost3D const& ghost_;
	mutable vector<Slope3D> rslopes_;
	mutable vector<Slope3D> naive_rslopes_;
	const bool slf_;
	const double shockratio_;
	const double diffusecoeff_;
	const double pressure_ratio_;
	const bool SR_;
	const vector<string> calc_tracers_;
	const string skip_key_;
	mutable vector<size_t> to_skip_;
	const bool pressure_calc_;

	LinearGauss3D(const LinearGauss3D& origin);

	LinearGauss3D& operator=(const LinearGauss3D& origin);
};
#endif //LINEAR_GAUSS3D_HPP