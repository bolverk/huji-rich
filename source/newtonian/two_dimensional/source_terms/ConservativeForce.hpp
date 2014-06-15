/*! \file ConservativeForce.hpp
  \brief Abstract class for conservative force's acceleration
  \author Elad Steinberg
*/

#ifndef CONSFORCE_HPP
#define CONSFORCE_HPP 1

#include "../SourceTerm.hpp"

//! \brief Physical acceleration
class Acceleration
{
public:
	/*!
	\brief Calculates the acceleration that the cell feels
	\param tess The tessellation
	\param cells The primitive cells
	\param point The index of the cell to calculate
	\param fluxes The vector of the fluxes
	\param point_velocity The velocities of the mesh points
	\param hbc The hydro boudnary conditions
	\param time The simulation time
	\param dt The time step
	\return The calculated acceleration
	*/
  virtual Vector2D Calculate(Tessellation const& tess,
			     vector<Primitive> const& cells,
			     int point,
			     vector<Conserved> const& fluxes,
			     vector<Vector2D> const& point_velocity,
			     HydroBoundaryConditions const& hbc,
			     vector<vector<double> > const& tracers,
			     double time,
			     double dt)=0;

  virtual ~Acceleration(void);
};
/*! \brief Class for conservative forces
\author Elad Steinberg
*/
class ConservativeForce: public SourceTerm
{
public:
	/*!
	\brief Class constructor
	\param acc The acceleration force
	\param DtCalc Should we calcualte the time step for the force
	*/
	ConservativeForce(Acceleration& acc,bool DtCalc=false);

	/*!
	\brief Class destructor
	*/
	~ConservativeForce(void);

	Conserved Calculate(Tessellation const& tess,
			    vector<Primitive> const& cells,
			    int point,
			    vector<Conserved> const& fluxes,
			    vector<Vector2D> const& point_velocity,
			    HydroBoundaryConditions const& hbc,
			    vector<vector<double> > const &tracer_extensive,
			    vector<double> &dtracer,vector<double> const& lengthes,
			    double time,double dt);

	/*!
	\brief Returns a the smallest time step for all cells based on dt=sqrt(R/g) where R is the cell's width and g is the acceleration
	\returns The time step
	*/
	double GetTimeStep(void) const;

private:
	Acceleration& acc_;
	bool DtCalc_;
	double dt_;
	double last_time_;

  ConservativeForce(const ConservativeForce& origin);
  ConservativeForce& operator=(const ConservativeForce& origin);
};

#endif // CONSFORCE_HPP
