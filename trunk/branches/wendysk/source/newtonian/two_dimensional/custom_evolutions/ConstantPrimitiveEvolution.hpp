#ifndef COSNTANTPRIMITIVEEVOLUTION_HPP
#define COSNTANTPRIMITIVEEVOLUTION_HPP 1

#include <fstream>
#include <string>
#include <iostream>
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
	\param mass_count A flag to whether track the mass influx into the cell
	\param n The number of cells with ConstantPrimitiveEvolution, the cells must be the first n in the cells vector
	*/
	ConstantPrimitiveEvolution(bool mass_count=false,int n=0);
	/*!
	\brief Class destructor
	*/
	~ConstantPrimitiveEvolution(void);

	Conserved CalcFlux(Tessellation const* tessellation,
		vector<Primitive> const& cells,	double dt,
		SpatialReconstruction* interpolation,Edge const& edge,
		Vector2D const& facevelocity,
		RiemannSolver const& rs,int index,
		HydroBoundaryConditions const* boundaryconditions,double time,
		vector<vector<double> > const& tracers); 

	Primitive UpdatePrimitive(vector<Conserved> const& conservedintensive,
		EquationOfState const* eos,vector<Primitive>& cells,int index);
	/*!
	\brief Returns the mass flux into the cell
	\return The mass flux
	*/
	double GetMassFlux(void)const;
	/*!
	\brief Reads the mass flux from a file
	\param loc The location of the file to read
	*/
	void ReadMassFlux(string loc)const;
	/*!
	\brief Writes the mass flux to a file
	\param loc The location of the file to write
	*/
	void WriteMassFlux(string loc)const;
private:
	bool mass_count_;
	double mass_flux,mass_fluxt;
	int N_;
};

#endif //COSNTANTPRIMITIVEEVOLUTION_HPP
