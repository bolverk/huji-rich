/*! \file ResetDump.hpp
  \brief Data needed to restart a simulation
  \author Elad Steinberg
 */

#ifndef RESTDUMP_HPP
#define RESTDUMP_HPP 1
#include "../common/hydrodynamic_variables.hpp"
#include <vector>

using std::vector;

//! \brief Container for grid and hydrodynamical cells
class HydroSnapshot
{
public:
  /*! \brief Class constructor
    \param mesh_points The mesh points
    \param cells The primitive cells
  */
  HydroSnapshot(vector<Vector2D> const& mesh_points,
		vector<Primitive> const& cells);
  //! \brief Default constructor
  HydroSnapshot();
  //! \brief The mesh points
  vector<Vector2D> mesh_points;
  //! \brief The primitive cells
  vector<Primitive> cells;
};

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
	//! \brief The points of the processors
	vector<Vector2D> procmesh;
	//! \brief The indeces of the custom evolution
	vector<size_t> cevolve;
};

#endif //RESTDUMP_HPP
