/*! \file CenterGravity.hpp
\brief Point source gravity force
\author Elad Steinberg
*/

#ifndef SIMPLEGRAVITY_HPP
#define SIMPLEGRAVITY_HPP 1

#include "../SourceTerm.hpp"

//! \brief Gravitational acceleration due to a pointlike mass
class SimpleGravity : public SourceTerm
{
public:
	/*!
	\brief Class constructor
	\param M The mass of the point source
	\param Rmin The softenning length
	\param center The location of the point source
	\param SoftStart The maximum radius at which the softening length works.
	*/
	SimpleGravity(double M, double Rmin, double SoftStart = 0, Vector2D center = Vector2D(),double omega=0);

	Conserved Calculate
		(Tessellation const& tess,
		const PhysicalGeometry& pg,
		vector<Primitive> const& cells,
		int point,
		vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const& hbc,
		vector<vector<double> > const& tracers,
		vector<double>& dtracer, vector<double> const& lengthes,
		double t,
		double dt);
	
	/*!
	\brief Returns a the smallest time step for all cells based on dt=sqrt(R/g) where R is the cell's width and g is the acceleration
	\returns The time step
	*/
	double GetTimeStep(void) const;
private:
	double M_;
	double Rmin_;
	double  softlength_;
	Vector2D _center;
	double omega_;
	double dt_;
	bool first_time_;
	double time_;
};

#endif // SIMPLEGRAVITY_HPP
