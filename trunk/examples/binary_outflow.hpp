#ifndef BINARY_OUTFLOW_HPP
#define BINARY_OUTFLOW_HPP 1

#include "../source/newtonian/two_dimensional/SourceTerm.hpp"
#include "../source/newtonian/common/equation_of_state.hpp"

/*!
	\brief Example of creating a custom force for Binary outflow/inflow
*/

class BinaryOutflow: public SourceTerm
{
private:
	double mdot_,e_,cs_;
	int nx_,ny_;
	int cell_index_;

public:
	/*!
	\brief Class constructor
	\param mdot The mass per unit time that is ejected
	\param e The thermal energy per unit mass of the ejecta
	\param cs The velocity of the ejecta
	\param nx The number of cells in the x direction
	\param ny The number of cells in the y direction
	*/
	BinaryOutflow(double mdot,double e,double cs,int nx,int ny);

	/*!
		\brief Initializes the outflow parameters
		\param tess The tessellation
		\param xc The x coordinate of the source
		\param yc The x coordinate of the source
	*/
	void Init(Tessellation const& tess,double xc,double yc);

	/*!
		\brief Returns the index of the outflow cell
		\returns The cell index
	*/
	int GetCellIndex(void);

	Conserved Calculate(Tessellation const* tess,
	  vector<Primitive> const& cells,int point,vector<Conserved> const& fluxes,
	  vector<Vector2D> const& point_velocity,HydroBoundaryConditions const*hbc,
	  vector<vector<double> > const& tracers,vector<double> &dtracer,double t,
	  double dt);
};

#endif // BINARY_OUTFLOW_HPP