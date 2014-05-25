#ifndef RESET_HPP
#define RESET_HPP 1
#include "hdsim2d.hpp"
#include "hydrodynamics_2d.hpp"
#include <fstream>
#include <string>

//! \brief Contains data required to restart a simulation run
class ResetDump
{
public:
	//! \brief Class constructor
	ResetDump();
	//! \brief Class destructor
	~ResetDump();
	//! \brief Clears the data from the class
	void clear(void);
	//! \brief The hydro variables and mesh points
	HydroSnapshot snapshot;
	//! \brief The tracers
	vector<vector<double> > tracers;
	//! \brief The simulation time
	double time;
	//! \brief The Courant number
	double cfl;
	//! \brief The simulation time step number
	int cycle;
	//! \brief Coldflows flag
	bool coldflows;
	//! \brief Densityfloor flag
	bool densityfloor;
	//! \brief Coldflows kinetic energy ratio
	double a;
	//! \brief Coldflows potential energy ratio
	double b;
	//! \brief The density for density floor
	double densitymin;
	//! \brief The pressure for density floor
	double pressuremin;
};
/*!
\brief Outputs the simulation data into a restart file
\param location The filename of the dump
\param sim The sim data
*/
void ResetOutput(string location,hdsim const& sim);
/*!
\brief Reads the simulation data from a restart file
\param location The filename of the dump
\param dump The dump class that is overwritten
\param eos The equation of state
*/
void ResetRead(string location,ResetDump &dump,EquationOfState const* eos);

#endif // RESET_HPP
