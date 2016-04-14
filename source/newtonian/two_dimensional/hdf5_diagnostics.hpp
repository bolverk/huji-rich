/*! \file hdf5_diagnostics.hpp
  \brief Simulation output to hdf5 file format
  \author Elad Steinberg
 */

#ifndef HDF5_DIAG
#define HDF5_DIAG 1

#include <H5Cpp.h>
#include <string>
#include "hdsim2d.hpp"
#include "../../misc/int2str.hpp"
#include "diagnostics.hpp"
#include "../../tessellation/Delaunay.hpp"

//! \brief Container for snapshot data
class Snapshot
{
public:

  //! \brief Default constructor
  Snapshot(void);

  //! \brief Copy constructor
  //! \param source Source
  Snapshot(const Snapshot& source);

  //! \brief Mesh points
  vector<Vector2D> mesh_points;

  //! \brief Computational cells
  vector<ComputationalCell> cells;

  //! \brief Time
  double time;

  //! \brief Cycle number
  int cycle;

#ifdef RICH_MPI
  //! \brief Locations of cpus
  vector<Vector2D> proc_points;
#endif

  //! \brief THe names of the tracers and stickers
  TracerStickerNames tracerstickernames;
};

/*! \brief Load snapshot data into memory
  \param fname File name
  \param mpioverride Flag for not reading mpi data when MPI is on
  \return Snapshot data
 */
Snapshot read_hdf5_snapshot
(const string& fname,bool mpioverride=false);

//! \brief Addition data to be written in a snapshot
class DiagnosticAppendix
{
public:

  /*! \brief Calculates additional data
    \param sim Hydrodynamic simulation
    \return Calculated data
   */
  virtual vector<double> operator()(const hdsim& sim) const = 0;

  /*! \brief Returns the name of the new field
   */
  virtual string getName(void) const = 0;

  //! \brief Class destructor
  virtual ~DiagnosticAppendix(void);
};

/*!
\brief Writes the simulation data into an HDF5 file
\param sim The hdsim class of the simulation
\param fname The name of the output file
\param appendices Additional data to be written to snapshot
*/
void write_snapshot_to_hdf5(hdsim const& sim,string const& fname,
			    const vector<DiagnosticAppendix*>& appendices=vector<DiagnosticAppendix*>());
/*!
\brief Reads an HDF5 snapshot file in order to restart the simulation
\param dump The dump data structure, should be when passed
\param fname The path to the HDF5 file
\param eos The equation of state
*/
void read_hdf5_snapshot(ResetDump &dump,string const& fname,EquationOfState
	const* eos);
/*!
\brief Converts an HDF5 snapshot file to the RICH custom reset binary format
\param input The path to the HDF5 file
\param output The path to the new binary file
*/
void ConvertHDF5toBinary(string const& input, string const& output);

/*!
\brief Writes the Delaunay triangulation data into an HDF5 file
\param tri The triangulation
\param filename The name of the output file
*/
void WriteDelaunay(Delaunay const& tri, string const& filename);

/*!
\brief Reads an HDF5 snapshot file in order to restart the simulation with a different cpu number
\return dump The snapshot structure relevent for current cpu
\param filename File name
\param proctess Tessellation of the processors
\param snapshot_number Number of old cpus
*/
Snapshot ReDistributeData(string const& filename, Tessellation const& proctess, size_t snapshot_number);

/*!
\brief Reads an HDF5 snapshot file in order to restart the simulation with a different cpu number. Slower but more robust then other method
\return dump The snapshot structure relevent for current cpu
\param filename File name
\param proctess Tessellation of the processors
\param snapshot_number Number of old cpus
\param mpioverride Flag for not reading mpi data when MPI is on
*/
Snapshot ReDistributeData2(string const& filename, Tessellation const& proctess, size_t snapshot_number,
	bool mpioverride = false);


/*!
\brief Writes the tessellation data into an HDF5 file
\param tess The tessellation
\param filename The name of the output file
*/
void WriteTess(Tessellation const& tess, string const& filename);
#endif // HDF5_DIAG
