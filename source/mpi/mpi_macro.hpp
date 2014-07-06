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
  \param second Second iterator
  \param val Value
  \return Iterator of the cell that contains val
  \todo Move this function to source/misc/utils.hpp
 */
template<class InputIterator, class T>
  InputIterator Find (InputIterator first, InputIterator last, const T& val)
{
  while (first!=last) {
    if (*first==val) return first;
    ++first;
  }
  return last;
}

// npoints includes ghost points
void PeriodicUpdateCells(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
	int npoints);

void PeriodicGradExchange(vector<ReducedPrimitiveGradient2D> &grad,
	vector<vector<int> > const& sentcells,int npoints);

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

int MPI_Send_Vector2D(Vector2D const& vec,int dest,int tag, MPI_Comm comm);

int MPI_VectorSend_Vector2D(vector<Vector2D> const& vec,int dest,int tag, MPI_Comm comm);

int MPI_Recv_Vector2D(Vector2D &vec,int source, int tag, MPI_Comm comm);

int MPI_VectorRecv_Vector2D(vector<Vector2D> &vec,int source, int tag, MPI_Comm comm);

int MPI_VectorBcast_Vector2D(vector<Vector2D> &vec,int root, MPI_Comm comm,int rank);

bool PointInsideCell(Tessellation const& tess,int cell_index,Vector2D const & point);

//void BuildTree(ANNkd_tree *&tree,ANNpointArray &treePoints,Tessellation const& tess);

//void KillTree(ANNkd_tree *&tree,ANNpointArray &treePoints);

//int FindContainingCell(ANNkd_tree *tree,Vector2D const& point);

void ConvertVector2DToDouble(vector<Vector2D> const& vec,vector<double> &res);

vector<Vector2D> MPI_MassSendRecvVectorVector2D(vector<vector<Vector2D> > const& tosend,
	vector<int> const& proclist,vector<int> const& procorder);

vector<int> GetProcOrder(int rank,int worldsize);

void SendRecvHydro(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> > sentcells,
	vector<int> sentprocs,EquationOfState const& eos,int totalpoints);

void SendRecvExtensive(vector<Conserved> const& cons,vector<vector<double> > const& 
	tracers,vector<size_t> const& customevolutions,vector<vector<int> > const& sentcells,
	vector<int> const& sentprocs,vector<Conserved> &ptoadd,vector<vector<double> >
	 &ttoadd,vector<size_t> &ctoadd);

void SendRecvOldVector2D(vector<Vector2D> const& points,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<Vector2D> &toadd);

int MPI_SendVectorPrimitive(vector<Primitive> const& vec,int dest,int tag,
	MPI_Comm comm);

void MPI_RecvVectorPrimitive(vector<Primitive> &vec,int dest,int tag,
	MPI_Comm comm,EquationOfState const& eos);

int MPI_SendVectorConserved(vector<Conserved> const& vec,int dest,int tag,
	MPI_Comm comm);

void MPI_RecvVectorConserved(vector<Conserved> &vec,int dest,int tag,
	MPI_Comm comm);


int MPI_SendVectorTracer(vector<vector<double> > const& vec,int dest,int tag,
	MPI_Comm comm);

void MPI_RecvVectorTracer(vector<vector<double> > &vec,int dest,int tag,
	MPI_Comm comm,int ntracer);

int MPI_SendVectorGrad(vector<ReducedPrimitiveGradient2D> const&vec,int dest,int
	tag,MPI_Comm comm);

void MPI_RecvVectorGrad(vector<ReducedPrimitiveGradient2D> &vec,int dest,int
	tag,MPI_Comm comm,int gradlength);

void SendRecvGrad(vector<ReducedPrimitiveGradient2D> &grads,
	vector<vector<int> > sentcells,vector<int> sentprocs,int totalpoints);

void SendRecvVelocity(vector<Vector2D> &vel,
	vector<vector<int> > sentcells,vector<int> sentprocs,int totalpoints);

void KeepLocalPoints(vector<Conserved> &cons,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<int> const& localpoints);

void SendRecvShockedStatus(vector<char> const& shockedcells,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<char> &btoadd);

void SendRecvVectorDouble(vector<double> const& vec,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<double> &toadd);

#endif

#endif //MPI_MACRO

