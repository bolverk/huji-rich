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


//! \brief Appendix to data dump
class DiagnosticAppendix3D
{
public:

	/*! \brief Calculates additional data
	\param sim Hydrodynamic simulation
	\return Calculated data
	*/
	virtual vector<double> operator()(const HDSim3D& sim) const = 0;

	/*! \brief Returns the name of the new field
	*/

	/*! \brief Get appendix title
	  \return Title
	 */
	virtual string getName(void) const = 0;

	//! \brief Class destructor
	virtual ~DiagnosticAppendix3D(void);
};

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

#ifdef RICH_MPI
	//! \brief Processors Mesh points
	vector<Vector3D> proc_points;
#endif

	//! \brief Volume of cells
	vector<double> volumes;

	//! \brief Computational cells
	vector<ComputationalCell3D> cells;

	//! \brief Time
	double time;

	//! \brief Cycle number
	int cycle;

	//! \brief THe names of the tracers and stickers
  pair<vector<string>, vector<string> > tracerstickernames;

	//! \brief Lower left corner
	Vector3D ll;

  //! \brief Upper right corner
	Vector3D ur;
};

#if RICH_MPI
/*! \brief Load snapshot data into memory
\param fname File name
\param mpi_write Flag for providing parallelisation data
\param fake_rank Process id
\return Snapshot data
*/
#else
/*! \brief Load snapshot data into memory
\param fname File name
\param fake_rank Process id
\return Snapshot data
*/
#endif // RICH_MPI
Snapshot3D ReadSnapshot3D(const string& fname
#ifdef RICH_MPI
	,bool mpi_write = false,int fake_rank=-1
#endif
);

#ifdef RICH_MPI
/*! \brief Redistribute data between the different processes
  \param filename Name of output file
  \param proctess Meta tessellation
  \param snapshot_number Number of snapshot
  \param mpi_write Parallel output flag
  \return Hydrodynamic snapshot
 */
Snapshot3D ReDistributeData3D(string const& filename, Tessellation3D const& proctess, size_t snapshot_number,bool mpi_write=false);
#endif


/*! \brief Write voronoi data to a file
  \param tri Voronoit tessellation
  \param filename Name of output file
 */
void WriteVoronoi(Voronoi3D const& tri, std::string const& filename);

std::vector<Vector3D> ReadVoronoiPoints(std::string const& filename);

#if RICH_MPI
/*! \brief Write snapshot to file
  \param sim Simulation
  \param filename name of output file
  \param appendices Custom fields
  \param mpi_write Determines whether to write parallisation data
 */
#else
/*! \brief Write snapshot to file
  \param sim Simulation
  \param filename name of output file
  \param appendices Custom fields
 */
#endif // RICH_MPI
void WriteSnapshot3D(HDSim3D const& sim, std::string const& filename,
	const vector<DiagnosticAppendix3D*>& appendices = vector<DiagnosticAppendix3D*>()
#ifdef RICH_MPI
	,bool mpi_write = true
#endif
);
#endif // HDF_WRITE
