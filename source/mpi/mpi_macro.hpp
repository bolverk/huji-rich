/*! \file mpi_macro.hpp
  \brief A collection of useful utilities for mpi
  \author Elad Steinberg
 */

#ifndef MPI_MACRO
#define MPI_MACRO 1

#define _USE_MATH_DEFINES 
#include <cmath>
#include "../newtonian/two_dimensional/OuterBoundary.hpp"
#include "../tessellation/tessellation.hpp"
#include "../tessellation/geotests.hpp"
#include "../treecode/ANN.h"
#include "../newtonian/common/hydrodynamic_variables.hpp"
#include "../newtonian/common/equation_of_state.hpp"
#include "../misc/utils.hpp"
#include "../newtonian/two_dimensional/ReducedPrimitiveGradient2D.hpp"

using namespace std;

/*! \brief Find the iterator of a certain value in an array
  \param first First iterator
  \param last Iterator pointing to last position
  \param val Value
  \return Iterator of the cell that contains val
  \todo Move this function to source/misc/utils.hpp
 */
template<class InputIterator, class T>
  InputIterator Find (InputIterator first,InputIterator last, const T& val)
{
  while (first!=last) {
    if (*first==val) return first;
    ++first;
  }
  return last;
}

  /*!
  \brief Creates the simulation data for the ghost points with periodic boundaries
  \param cells The primitive cells
  \param tracers The intensive tracers
  \param customevolutions The indeces of the custom evolutions
  \param sentcells The indeces of the cells that are ghost points (ordered with the same order as they were added to the tessellation)
  \param npoints The total number of mesh points (including the ghost cells)
  */
// npoints includes ghost points
void PeriodicUpdateCells(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
	int npoints);
/*!
  \brief Copies the relevant gradiant to the ghost points with periodic boundaries
  \param grad The original gradients
  \param sentcells The indeces of the cells that are ghost points (ordered with the same order as they were added to the tessellation)
  \param npoints The total number of mesh points (including the ghost cells)
  */
void PeriodicGradExchange(vector<ReducedPrimitiveGradient2D> &grad,
	vector<vector<int> > const& sentcells,int npoints);
/*!
  \brief Copies the relevant velocity of the mesh points to the ghost points with periodic boundaries
  \param vel The original velocities
  \param sentcells The indeces of the cells that are ghost points (ordered with the same order as they were added to the tessellation)
  \param npoints The total number of mesh points (including the ghost cells)
  */
void PeriodicVelocityExchange(vector<Vector2D> &vel,
	vector<vector<int> > const& sentcells,int npoints);

#ifdef RICH_MPI

#include <mpi.h>

/*! \brief Returns the mpi rank of the current process
  \return MPI rank
 */
int get_mpi_rank(void);

/*! \brief Returns the mpi size of the current group
  \return MPI size
 */
int get_mpi_size(void);
/*!
\brief MPI wrapper to send a single Vector2D
\param vec The Vector2D to send
\param dest The id of the receiving cpu
\param tag The tag of the msg
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_Send_Vector2D(Vector2D const& vec,int dest,int tag, MPI_Comm comm);
/*!
\brief MPI wrapper to send a vector of Vector2D
\param vec The vector of Vector2D to send
\param dest The id of the receiving cpu
\param tag The tag of the msg
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_VectorSend_Vector2D(vector<Vector2D> const& vec,int dest,int tag, MPI_Comm comm);
/*!
\brief MPI wrapper to receive a single Vector2D
\param vec The Vector2D that will receive the data
\param source The id of the sending cpu
\param tag The tag of the msg
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_Recv_Vector2D(Vector2D &vec,int source, int tag, MPI_Comm comm);
/*!
\brief MPI wrapper to receive a vector of Vector2D
\param vec The vector of Vector2D that will receive the data
\param source The id of the sending cpu
\param tag The tag of the msg
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_VectorRecv_Vector2D(vector<Vector2D> &vec,int source, int tag, MPI_Comm comm);
/*!
\brief MPI wrapper to broadcast a vector of Vector2D
\param vec The vector Vector2D that will receive/send the data
\param root The id of the broadcasting cpu
\param comm The MPI communicator
\param rank The rank of the current cpu
\return The MPI error code
*/
int MPI_VectorBcast_Vector2D(vector<Vector2D> &vec,int root, MPI_Comm comm,int rank);
/*!
\brief Checks if a point is inside a Voronoi cell
\param tess The tessellation
\param cell_index The index of the cell to check
\param point The point to cehck
\return True of false
*/
bool PointInsideCell(Tessellation const& tess,int cell_index,Vector2D const & point);

//void BuildTree(ANNkd_tree *&tree,ANNpointArray &treePoints,Tessellation const& tess);

//void KillTree(ANNkd_tree *&tree,ANNpointArray &treePoints);

//int FindContainingCell(ANNkd_tree *tree,Vector2D const& point);

/*!
\brief Converts a vector of Vector2D into a vector of doubles
\param vec The vector of Vector2D
\param res The vector of doubles given as output
*/
void ConvertVector2DToDouble(vector<Vector2D> const& vec,vector<double> &res);

/*!
\brief Sends and receives data from several cpus, each cpu sends/recv a vector of Vector2D
\param tosend The data to send, outermost vector labels the cpu
\param proclist The id of cpus to talk with
\param procorder The order in which to talk with processors
\return The received data
*/
vector<Vector2D> MPI_MassSendRecvVectorVector2D(vector<vector<Vector2D> > const& tosend,
	vector<int> const& proclist,vector<int> const& procorder);

/*!
\brief Returns the communication order
\param rank The id of the current cpu
\param worldsize The number of total cpus
\return The communication order
*/
vector<int> GetProcOrder(int rank,int worldsize);

/*!
\brief Sends and receives the hydro for ghost cells
\param cells The primitive cells
\param tracers The intensive tracers
\param customevolutions The custom evolution indeces
\param sentcells The indeces of the cells to send
\param sentprocs The ids of the cpus to talk with
\param eos The equation of state
\param Nghost The indeces of each ghost point in the vector of points that the tessellation holds
\param totalpoints The total number of points in teh tessellation (including ghost)
*/
void SendRecvHydro(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> >const& sentcells,
	vector<int>const& sentprocs,EquationOfState const& eos,
	vector<vector<int> > const& Nghost,int totalpoints);
/*!
\brief Sends and receives the extensive hydro, typically used for excahgning real (not ghost) points that are transfered between cpus
\param cons The extensive cells
\param tracers The extensive tracers
\param customevolutions The custom evolution indeces
\param sentcells The indeces of the cells to send
\param sentprocs The ids of the cpus to talk with
\param ptoadd The received extensive cells
\param ttoadd The received extensive tracers
\param ctoadd The received custom evolution indeces
*/
void SendRecvExtensive(vector<Conserved> const& cons,vector<vector<double> > const& 
	tracers,vector<size_t> const& customevolutions,vector<vector<int> > const& sentcells,
	vector<int> const& sentprocs,vector<Conserved> &ptoadd,vector<vector<double> >
	 &ttoadd,vector<size_t> &ctoadd);
/*!
\brief Send/recv the old mesh generating points
\param points The mesh generating points
\param sentcells The indeces of the cells to send
\param sentprocs Process id's to which the message should be sent
\param toadd The received mesh generating points
*/
void SendRecvOldVector2D(vector<Vector2D> const& points,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<Vector2D> &toadd);

/*!
\brief send a vector of primitives
\param vec The vector to send
\param dest The id of the cpu to recv
\param tag The MPI msg tag
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_SendVectorPrimitive(vector<Primitive> const& vec,int dest,int tag,
	MPI_Comm comm);
/*!
\brief recv a vector of primitives
\param vec The vector that recveives
\param dest The id of the cpu taht sent the data
\param tag The MPI msg tag
\param comm The MPI communicator
\param eos The equation of state
*/
void MPI_RecvVectorPrimitive(vector<Primitive> &vec,int dest,int tag,
	MPI_Comm comm,EquationOfState const& eos);
/*!
\brief send a vector of conserved variables
\param vec The vector to send
\param dest The id of the cpu to recv
\param tag The MPI msg tag
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_SendVectorConserved(vector<Conserved> const& vec,int dest,int tag,
	MPI_Comm comm);
/*!
\brief recv a vector of conserved
\param vec The vector that recveives
\param dest The id of the cpu taht sent the data
\param tag The MPI msg tag
\param comm The MPI communicator
*/
void MPI_RecvVectorConserved(vector<Conserved> &vec,int dest,int tag,
	MPI_Comm comm);
/*!
\brief send a vector of tracers
\param vec The vector to send
\param dest The id of the cpu to recv
\param tag The MPI msg tag
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_SendVectorTracer(vector<vector<double> > const& vec,int dest,int tag,
	MPI_Comm comm);
/*!
\brief recv a vector of tracers
\param vec The vector that recveives
\param dest The id of the cpu taht sent the data
\param tag The MPI msg tag
\param comm The MPI communicator
\param ntracer The number of tracers
*/
void MPI_RecvVectorTracer(vector<vector<double> > &vec,int dest,int tag,
	MPI_Comm comm,int ntracer);
/*!
\brief Send a vector of gradients
\param vec The vector to send
\param dest The id of the cpu taht sent the data
\param tag The MPI msg tag
\param comm The MPI communicator
\return The MPI error code
*/
int MPI_SendVectorGrad(vector<ReducedPrimitiveGradient2D> const&vec,int dest,int
	tag,MPI_Comm comm);
/*!
\brief Recv a vector of gradients
\param vec The vector that receives
\param dest The id of the cpu taht sent the data
\param tag The MPI msg tag
\param comm The MPI communicator
\param gradlength How many doubles does the gradient have? (Used for determining tracer number)
\return The MPI error code
*/
void MPI_RecvVectorGrad(vector<ReducedPrimitiveGradient2D> &vec,int dest,int
	tag,MPI_Comm comm,int gradlength);
/*!
\brief Send/Recv gradients of ghost cells
\param grads The original vector of gradients
\param sentcells The indeces of the gradients to send
\param sentprocs The ids of the cpus to talk with
\param Nghost The indeces of each ghost point in the vector of points that the tessellation holds
\param totalpoints The total number of points in the tessellation (including ghost)
*/
void SendRecvGrad(vector<ReducedPrimitiveGradient2D> &grads,
	vector<vector<int> >const& sentcells,vector<int> sentprocs,
	vector<vector<int> > const& Nghost,int totalpoints);
/*!
\brief Send/Recv velocities of ghost cells
\param vel The original vector of velocities
\param sentcells The indeces of the gradients to send
\param sentprocs The ids of the cpus to talk with
\param Nghost The indeces of each ghost point in the vector of points that the tessellation holds
\param totalpoints The total number of points in the tessellation (including ghost)
*/
void SendRecvVelocity(vector<Vector2D> &vel,
	vector<vector<int> >const& sentcells,vector<int> sentprocs,
	vector<vector<int> > const& Nghost,int totalpoints);
/*!
\brief Resizes and modifies the vectors to retain only cells that are inside the local cpu domain
\param cons The extensive cells
\param tracers The extensive tracers
\param customevolutions The custom evolution indeces
\param localpoints The indeces of the points to keep
*/
void KeepLocalPoints(vector<Conserved> &cons,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<int> const& localpoints);
/*!
\brief Send/Recv the shock state of cells
\param shockedcells The shocked state of the cells
\param sentcells The indeces of the vector to send
\param sentprocs The ids of the cpus to talk with
\param btoadd The received result
*/
void SendRecvShockedStatus(vector<char> const& shockedcells,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<char> &btoadd);
/*!
\brief Send/Recv a vector of doubles
\param vec the vector to send
\param sentcells The indeces of the vector to send
\param sentprocs The ids of the cpus to talk with
\param toadd The received result
*/
void SendRecvVectorDouble(vector<double> const& vec,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<double> &toadd);

/*! 
\brief Send/Recv the list of points to remove that are on the edges between cpus
\param GhostIndeces The indeces of the boundary points that were sent to current rank (this is recieved)
\param BoundaryPoints The points in the boundary that current cpu sent, each vector should be sorted
\param SentPoints The list of points that were sent to neighboring cpus
\param SentProcs The list of neighboring cpus
*/
void SendRecvGhostIndeces(vector<vector<int> > &GhostIndeces,vector<int>
	const& BoundaryPoints,vector<vector<int> > const& SentPoints,vector<int> const&
	SentProcs);

#endif

#endif //MPI_MACRO

