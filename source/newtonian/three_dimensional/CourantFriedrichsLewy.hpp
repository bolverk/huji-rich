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
	  \param no_calc Name of stickers of cells not to take into account for dt calculation
	 */
	explicit CourantFriedrichsLewy(double cfl, double SourceCFL, SourceTerm3D const& source,
		std::vector<std::string> no_calc = std::vector<std::string> (),	bool debug = false);

	double operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, const EquationOfState& eos,
		const vector<Vector3D>& face_velocities, const double time) const;

	void SetTimeStep(double dt);

private:
	const double cfl_, sourcecfl_;
	SourceTerm3D const& source_;
	std::vector<std::string> const no_calc_;
	const bool debug_;
	mutable bool first_try_;
	mutable double dt_first_;
	mutable double last_time_;
};

#endif // COURANT_FRIEDRICHS_LEWY_HPP
