/*! \file ConstantPrimitiveEvolution.hpp
  \brief Keeps primitive variables in a cell constant
  \author Elad Steinberg
  \todo Change file name to underscore seperated lowercase words, i.e. constant_primitive_evolution.hpp
 */

#ifndef COSNTANTPRIMITIVEEVOLUTION_HPP
#define COSNTANTPRIMITIVEEVOLUTION_HPP 1

#include <string>
#include "../CustomEvolution.hpp"
#include "../hydrodynamics_2d.hpp"

/*!
\brief Custom evolution for cells that keeps the cell constant.
\author Elad Steinberg
*/
class ConstantPrimitiveEvolution: public CustomEvolution
{
public:
	/*!
	\brief Class constructor
	\param prim The primitive to keep constant
	\param tracer The tracer to keep constant
	\param mass_count A flag to whether track the mass influx into the cell
	\param n The number of cells with ConstantPrimitiveEvolution, the cells must be the first n in the cells vector
	*/
	ConstantPrimitiveEvolution(Primitive const& prim,vector<double> const& tracer,
		bool mass_count=false,int n=0);
	/*!
	\brief Class destructor
	*/
	~ConstantPrimitiveEvolution(void);

	Conserved CalcFlux(Tessellation const& tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction& interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const& boundaryconditions,double time,
		vector<vector<double> > const& tracers);

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const& eos,vector<Primitive>& cells,int index,
		Tessellation const& tess,double time,vector<vector<double> > const& tracers);

	vector<double> UpdateTracer(int index,vector<vector<double> >
		const& tracers,vector<vector<double> > const& /*tracerchange*/,vector<Primitive> const& cells,Tessellation const& tess,
		double time);

	vector<double> CalcTracerFlux(Tessellation const& tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double dm,Edge const& edge,int index,double dt,double time,
		SpatialReconstruction const& interp,Vector2D const& vface);

	bool ShouldForceTracerReset(void)const;

	bool TimeStepRelevant(void)const;

	/*!
	\brief Returns the mass flux into the cell
	\return The mass flux
	*/
	double GetMassFlux(void)const;
	/*!
	\brief Sets the mass flux from a file
	\param m The total mass flux to date
	*/
	void SetMassFlux(double m);

	bool isRelevantToInterpolation(void) const;

private:
	Primitive prim_;
	vector<double> tracer_;
	bool mass_count_;
	double mass_flux,mass_fluxt;
	int N_;
};

#endif //COSNTANTPRIMITIVEEVOLUTION_HPP