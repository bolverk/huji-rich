/*! \file CourantFriedrichsLewy.hpp
  \brief Calculates the time step according to the CFL criterion
  \author Almog Yalinewich
 */

#ifndef COURANT_FRIEDRICHS_LEWY_HPP
#define COURANT_FRIEDRICHS_LEWY_HPP 1

#include "time_step_function3D.hpp"
#include "SourceTerm3D.hpp"

 //! \brief Calculates the time step according to the CFL criterion
class CourantFriedrichsLewy : public TimeStepFunction3D
{
public:

	/*! \brief Class constructor
	  \param cfl CFL number
	  \param source The source term for the simulation
	  \param SourceCFL CFL number for the source term
	 */
	explicit CourantFriedrichsLewy(double cfl,double SourceCFL, SourceTerm3D const& source,
		bool debug=false);

	double operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, const EquationOfState& eos,
		const vector<Vector3D>& face_velocities,const double time, TracerStickerNames const& tracerstickernames) const;

	void SetTimeStep(double dt);

private:
	const double cfl_,sourcecfl_;
	SourceTerm3D const& source_;
	const bool debug_;
	mutable bool first_try_;
	double dt_first_;
};

#endif // COURANT_FRIEDRICHS_LEWY_HPP
