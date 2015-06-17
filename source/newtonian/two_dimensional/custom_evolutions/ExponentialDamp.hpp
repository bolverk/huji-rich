/*! \file ExponentialDamp.hpp
\brief Damps primitve exponentialy to initial conditions
\author Elad Steinberg
*/

#ifndef ExponentialDamp_HPP
#define ExponentialDamp_HPP 1

#include <string>
#include "../CustomEvolution.hpp"
#include "../hydrodynamics_2d.hpp"

class ExponentialDamp : public CustomEvolution
{
public:
	/*!
	\brief Class constructor
	*/
	ExponentialDamp(double Rmin, double Rmax, double tau, SpatialDistribution const& density,
		SpatialDistribution const& pressure,SpatialDistribution const& xvel,SpatialDistribution const& yvel,
		EquationOfState const& eos);
	/*!
	\brief Class destructor
	*/
	~ExponentialDamp(void);

	Conserved CalcFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells, double dt,
		SpatialReconstruction& interpolation, Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs, int index,
		HydroBoundaryConditions const& boundaryconditions, double time,
		vector<vector<double> > const& tracers);

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const& eos, vector<Primitive>& cells, int index,
		Tessellation const& tess, double time, vector<vector<double> > const& tracers);

	vector<double> UpdateTracer(int index, vector<vector<double> >
		const& tracers, vector<vector<double> > const& /*tracerchange*/, vector<Primitive> const& cells, Tessellation const& tess,
		double time);

	vector<double> CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells, vector<vector<double> > const& tracers,
		double dm, Edge const& edge, int index, double dt, double time,
		SpatialReconstruction const& interp, Vector2D const& vface);

	bool ShouldForceTracerReset(void)const;

	bool TimeStepRelevant(void)const;

	bool isRelevantToInterpolation(void) const;

private:
	const double a_,b_,tau_;
	SpatialDistribution const& density_;
	SpatialDistribution const& pressure_;
	SpatialDistribution const& xvel_;
	SpatialDistribution const& yvel_;
	EquationOfState const& eos_;
};

#endif //ExponentialDamp_HPP
