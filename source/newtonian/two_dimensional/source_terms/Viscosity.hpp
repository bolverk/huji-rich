/*! \file AlphaDisc.hpp
\brief Alpha model for disc viscosity, assumes center of disc is at origin
\author Elad Steinberg
*/

#ifndef VISCOSITY_HPP
#define VISCOSITY_HPP 1

#include "../SourceTerm.hpp"

//! \brief Gravitational acceleration due to a pointlike mass
class Viscosity : public SourceTerm
{
public:
	/*!
	\brief Class constructor
	\param M The mass of the point source
	\param Rmin The softenning length
	\param center The location of the point source
	\param SoftStart The maximum radius at which the softening length works.
	*/
	Viscosity(double nu, SpatialReconstruction &grad);

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

	virtual double GetNu(Tessellation const& tess, const PhysicalGeometry& pg, vector<Primitive> const& cells, int point)const;

	double nu_;
	SpatialReconstruction &grads_;
};

#endif // VISCOSITY_HPP
