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

	//! \brief Volume of cells
	vector<double> volumes;

	//! \brief Computational cells
	vector<ComputationalCell3D> cells;

	//! \brief Time
	double time;

	//! \brief Cycle number
	int cycle;

#ifdef RICH_MPI
	//! \brief Locations of cpus
	vector<Vector3D> proc_points;
#endif

	//! \brief THe names of the tracers and stickers
	TracerStickerNames tracerstickernames;
};

/*! \brief Load snapshot data into memory
\param fname File name
\param mpioverride Flag for not reading mpi data when MPI is on
\return Snapshot data
*/
Snapshot3D ReadSnapshot3D(const string& fname, bool mpioverride = false);

void WriteVoronoi(Voronoi3D const& tri, std::string const& filename);

void WriteSnapshot3D(HDSim3D const& sim, std::string const& filename);
#endif // HDF_WRITE
