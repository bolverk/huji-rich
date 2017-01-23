/*! \file hdf_write.hpp
\brief Simulation output to hdf5 file format
\author Elad Steinberg
*/

#ifndef HDF_WRITE
#define HDF_WRITE 1

#include <H5Cpp.h>
#include <string>
#include "Voronoi3D.hpp"
#include "../../newtonian/three_dimensional/hdsim_3d.hpp"

//! \brief Container for snapshot data
class Snapshot3D
{
public:

	//! \brief Default constructor
	Snapshot3D(void);

	//! \brief Copy constructor
	//! \param source Source
	Snapshot3D(const Snapshot3D& source);

	//! \brief Mesh points
	vector<Vector3D> mesh_points;
};

void WriteVoronoi(Voronoi3D const& tri, std::string const& filename);

void WriteSnapshot(HDSim3D const& sim, std::string const& filename);
#endif // HDF_WRITE
