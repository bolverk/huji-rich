/*! \file static_voronoi_mesh.hpp
  \brief A static voronoi tessellation
  \author Almog Yalinewich
 */

#ifndef STATIC_VORONOI_MESH_HPP
#define STATIC_VORONOI_MESH_HPP 1

#define _USE_MATH_DEFINES
#include <cmath>
#include "Delaunay.hpp"
#include <list>
#include "tessellation.hpp"
#include "../misc/utils.hpp"
#include "../newtonian/two_dimensional/HydroBoundaryConditions.hpp"
#include "../newtonian/two_dimensional/RefineStrategy.hpp"
#include "../newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "voronoi_logger.hpp"
#include "ConvexHull.hpp"
#include "../misc/int2str.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#include "../mpi/mpi_macro.hpp"
#endif

using std::string;

//! \brief Voronoi tessellation class
class StaticVoronoiMesh : public Tessellation
{
public:
	/*!
	\brief Fixes the number in the NewIndex for periodic boundary conditions
	\param V The Voronoi
	\param other The real index of the cell whose edge we fix
	\param ToRefine The point that is refined
	\param NewIndex The index of the new point
	\param NewPoint The new point
	\return The index of the duplicated cor
	*/
	friend int FixPeriodNeighbor(StaticVoronoiMesh &V,int other,int ToRefine,
		int NewIndex,Vector2D const& NewPoint);

	Vector2D CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL,
		Vector2D rR,Vector2D f)const;

	void RemoveCells(vector<int> &ToRemove,vector<vector<int> > &VolIndex,
		vector<vector<double> > &Volratio);

	void RefineCells(vector<int> const& ToRefine,
		vector<Vector2D> const& directions,double alpha=1e-3);

	Tessellation* clone(void)const;

	vector<int> GetNeighbors(int index)const;

	/*!
	\brief Returns the list of neighbors including ghost points
	\return The neighbors
	\param index The cell whose neighbors are returned
	*/
	vector<int> GetLiteralNeighbors(int index)const;

	int GetOriginalIndex(int point) const;

#ifdef RICH_MPI
	void Initialise(vector<Vector2D> const& points,Tessellation const& vproc,
		OuterBoundary const* outer);
#endif

	void Initialise(vector<Vector2D> const& points,OuterBoundary const* bc);

	/*!
	\brief Class constructor
	\param points The points to build for
	\param bc The geometric boundary conditions
	*/
  StaticVoronoiMesh(vector<Vector2D> const& points,OuterBoundary const& bc);

	/*!
		\brief Copy constructor
		\param other The Voronoi to copy
	*/
  StaticVoronoiMesh(StaticVoronoiMesh const& other);

  #ifdef RICH_MPI
	void Update(vector<Vector2D> const& points,Tessellation const& vproc);
#endif // RICH_MPI

	void Update(vector<Vector2D> const& points);

  ~StaticVoronoiMesh(void);
	/*! \brief Get Total number of mesh generating points
    \return Number of mesh generating points
   */
  int GetPointNo(void) const;

  /*! \brief Returns Position of mesh generating point
    \param index Mesh generating point index
    \return Position of mesh generating point
   */
  Vector2D GetMeshPoint(int index) const;

  /*! \brief Returns the total number of faces
    \return Total number of faces
   */
  int GetTotalSidesNumber(void) const;

  /*! \brief Returns edge (interface between cells)
    \param index Face index
    \return Interface between cells
   */
  Edge const& GetEdge(int index) const;

  /*! \brief Returns the effective width of a cell
    \param index Cell index
    \return Effective cell width
   */
  double GetWidth(int index) const;

  /*! \brief Returns the volume of a cell
    \param index Cell index
    \return Cell volume
   */
  double GetVolume(int index) const;

//! \brief Returns a reference to a vector<int> containing the indexes of the edges related to a cell. \param index The index of the cell. \returns The reference to vector<int>
  vector<int>const& GetCellEdges(int index) const;

/*! \brief Returns Position of Cell's CM
    \param index Mesh generating point index (the cell's index)
    \return Position of CM
   */
  Vector2D const& GetCellCM(int index) const;

  vector<Vector2D>& GetMeshPoints(void);

  /*! \brief Outputs the grid data
  \param filename The path to the output file
  */
  void output(string filename);

  bool NearBoundary(int index) const;

  //! \brief Diagnostics method
  voronoi_loggers::VoronoiLogger* logger;

  /*! \brief Returns a reference to a list of all edges
    \return List of all edges
   */
  vector<Edge>& GetAllEdges(void);

  /*!
  \brief Calculates and send the ghost points that are needed for AMR with MPI
  \param ToRemove The local list of points to remove
  \param BoundaryRemove The list of points that are border points per proc that are removed
  \param BoundaryNeigh The indeces in the Nghost list for each neighbor of the points in BoundaryRemove
  */
  void FindBoundaryRemoveSend(vector<int> const& ToRemove,vector<vector<int> > &BoundaryRemove,
		vector<vector<vector<int> > > &BoundaryNeigh);

  vector<vector<int> >& GetDuplicatedPoints(void);

  vector<vector<int> >const& GetDuplicatedPoints(void)const;

  vector<int> GetDuplicatedProcs(void)const;

  vector<vector<int> >const& GetSentPoints(void)const;

  vector<int> GetSentProcs(void)const;

  vector<vector<int> >& GetGhostIndeces(void);

  vector<vector<int> >const& GetGhostIndeces(void)const;

  vector<int> GetSelfPoint(void)const;

  int GetTotalPointNumber(void)const;
  
  const vector<Edge>& getAllEdges(void) const;

  vector<Vector2D>& GetAllCM(void);

  void GetNeighborNeighbors(vector<int> &result, int point)const;
private:
  vector<Vector2D> vertices_;
	double eps;
	OuterBoundary const* obc;
	vector<Edge> cell_edges;
	vector<Edge> edges;
	vector<Vector2D> CM;
	vector<vector<int> > mesh_vertices; // Which edges does each mesh point have
	Delaunay Tri;
	vector<int> GhostProcs;
	vector<vector<int> > GhostPoints;
	vector<int> SentProcs;
	vector<vector<int> > SentPoints;
	vector<int> selfindex;
	vector<vector<int> > NGhostReceived;
	vector<vector<int> > OrgCorner;
	int Nextra;

	bool legal_edge(Edge *e);
	void build_v(void);//Builds the voronoi mesh
	StaticVoronoiMesh& operator=(const StaticVoronoiMesh& origin);
	Vector2D CalcCellCM(size_t index) const;
	void FindIntersectingPoints(vector<Edge> const& box_edges,
		vector<vector<int> > &toduplicate);
	vector<int> CellIntersectBoundary(vector<Edge> const&box_edges,int cell);
	void GetAdditionalBoundary(vector<vector<int> > &copied,
		vector<vector<int> > &result,vector<vector<int> > &totest);
	void GetCorners(vector<vector<int> > &copied,vector<vector<int> > &result);
	vector<int> AddPointsAlongEdge(size_t point,vector<vector<int> > const&copied,int side);
	void GetRealNeighbor(vector<int> &result,int point)const;
	vector<int> GetBorderingCells(vector<int> const& copied,
		vector<int> const& totest,int tocheck,vector<int> tempresult,int outer);
	bool CloseToBorder(int point,int &border);
	void GetToTest(vector<vector<int> > &copied,vector<vector<int> > &totest);
	void ConvexEdgeOrder(void);
	vector<int> FindEdgeStartConvex(int point);
	void SendRecv(vector<int> const& procorder,vector<int> const&
		proclist,vector<vector<int> > &data);
	void NonSendBoundary(vector<int> &	proclist,vector<vector<int> > &
		data,Tessellation const& v,vector<vector<int> > &totest,
		vector<Edge> const& boxedges);
	void NonSendCorners(vector<int> & proclist,vector<vector<int> > & data,
		Tessellation const& v);
	void RigidBoundaryPoints(vector<int> &points,Edge const& edge);
	void PeriodicBoundaryPoints(vector<int> &points,int edge_number);
	void CornerBoundaryPoints(vector<int> &points,int edge_number);
	boost::array<double,4> FindMaxCellEdges(void);
	vector<int> CellIntersectOuterBoundary(vector<Edge> const&box_edges,int cell);
	void FindIntersectingOuterPoints(vector<Edge> const&bedge,vector<vector<int> >
		&boxduplicate,vector<vector<int> > const&firstduplicated);
	void SendRecvRemove(vector<int> const& procorder,vector<int> const&
		proclist,vector<vector<int> > &data);
	void GetNeighborNeighborsMPI(vector<int> &result,int point);
};
/*! \brief Checks if a point is inside a Voronoi cell
\param cpoints The points of the cell in convex order
\param vec The point to check
\return True if the point is inside the cell, false otherwise
*/
bool PointInCell(vector<Vector2D> const& cpoints,Vector2D const& vec);

#endif // STATIC_VORONOI_MESH_HPP
