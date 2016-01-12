/*! \file Delaunay.hpp
  \brief Delaunay triangulation with mpi
  \author Elad Steinberg
*/

#ifndef DELAUNAY_HPP
#define DELAUNAY_HPP 1
#include "facet.hpp"
#include <algorithm>
#include <vector>
#include <stack>
#include <climits>
#include "geometry.hpp"
#include "../misc/universal_error.hpp"
#include "HilbertOrder.hpp"
#include "geotests.hpp"
#include "../misc/utils.hpp"
#include "delaunay_logger.hpp"
#include "Edge.hpp"
#include "../newtonian/two_dimensional/OuterBoundary.hpp"
#include "shape_2d.hpp"
#include "../misc/int2str.hpp"
#include "../newtonian/two_dimensional/diagnostics.hpp"
#ifdef RICH_MPI
#include "find_affected_cells.hpp"
#endif
/*! \brief The Delaunay data structure. Gets a set of points and constructs the Delaunay tessellation.
  \author Elad Steinberg
*/
class Delaunay
{
private:

#ifdef RICH_MPI
	vector<int> OrgIndex;

  int findSomeOuterPoint(void);

  pair<vector<vector<int> >, vector<vector<int> > > findOuterPoints(const Tessellation& t_proc,
	  const vector<Edge>& edge_list,const vector<Edge>& box_edges,vector<vector<int> > &NghostIndex);

  pair<vector<vector<int> >, vector<int> > FindOuterPoints2
  (const Tessellation& t_proc,
   const vector<Edge>& edge_list,
   vector<vector<int> > &to_duplicate,
   vector<vector<int> >& self_points,
   const vector<Edge>& box_edges,
	  vector<vector<int> > &NghostIndex);

  vector<vector<int> > boundary_intersection_check
  (const vector<Edge>& edges,
   const vector<vector<int> >& to_duplicate);
#endif // RICH_MPI

  enum Sides{RIGHT,UP,LEFT,DOWN,LU,LD,RU,RD};
  int lastFacet; //last facet to be checked in Walk
  bool CalcRadius;

  class DataOnlyForBuild
  {
  public:
    DataOnlyForBuild();
    DataOnlyForBuild(DataOnlyForBuild const& other);
    DataOnlyForBuild& operator=(DataOnlyForBuild const& other);
    vector<vector<char> > copied;
  };

  vector<double> radius;
  vector<Vector2D> cell_points;
  bool PointWasAdded;
  int last_facet_added;
  vector<facet> f;
  vector<Vector2D> cor;
  int length;
  size_t olength;
  int location_pointer;
  int last_loc;

  bool IsOuterFacet(int facet)const;
  void add_point(size_t index,stack<std::pair<size_t, size_t> > &flip_stack);
  void flip(size_t i,size_t j, stack<std::pair<size_t, size_t> > & flip_stack);
  size_t Walk(size_t point);
  void CheckInput();
  double CalculateRadius(int facet);
  int FindPointInFacet(int facet,int point);
  double FindMaxRadius(int point);
  void FindContainingTetras(int StartTetra,int point,vector<int> &tetras);
  vector<int> FindContainingTetras(int StartTetra, int point);
  vector<vector<int> > FindOuterPoints(vector<Edge> const& edges);
  vector<vector<int> > FindOuterPointsMPI(OuterBoundary const* obc,vector<Edge> const& edges,
	  Tessellation const& tproc,vector<vector<int> > &Nghost,vector<int> &proclist);
  bool IsTripleOut(int index) const;
  int FindTripleLoc(facet const& f)const;
  void AddFacetDuplicate(int index,vector<vector<int> > &toduplicate,vector<Edge>
	const& edges,vector<bool> &checked)const;
  void AddOuterFacets(int tri,vector<vector<int> > &toduplicate,vector<Edge>
	const& edges,vector<bool> &checked);

  vector<vector<int> > AddOuterFacetsMPI
  (int point,
   vector<vector<int> > &toduplicate,
   vector<int> &neigh,
   vector<bool> &checked,
   const Tessellation& tproc,
   const vector<Edge>& own_edges,
   bool recursive = false);

  void AddRigid(vector<Edge> const& edges,
	vector<vector<int> > &toduplicate);
  vector<vector<int> > AddPeriodic(OuterBoundary const* obc,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate);
  void AddHalfPeriodic(OuterBoundary const* obc,vector<Edge> const& edges,
	vector<vector<int> > &toduplicate);
  double GetMaxRadius(int point,int startfacet);
  void SendRecvFirstBatch(vector<vector<Vector2D> > &tosend,
	  vector<int> const& neigh,vector<vector<int> > &Nghost);
  vector<int> GetOuterFacets(int cur_facet,int real_point,int olength);

  Delaunay& operator=(const Delaunay& origin);

public:

#ifdef RICH_MPI
  /*! \brief Retrieves the original index of a point (in case a point was duplicated)
    \param index Index of a point
    \return Original index
   */
	int GetOrgIndex(int index)const;
#endif
  /*! \brief Changes the cor olength
    \param n The new length;
  */
  void ChangeOlength(int n);

  /*! \brief Changes the cor length
    \param n The new length
  */
  void Changelength(int n);

   /*! \brief Allows to change the cor
    \return Refrence to the cor vector
  */
  vector<Vector2D>& ChangeCor(void);

  /*! \brief Access to coordinates
    \return Reference to coordinates vector
   */
  const vector<Vector2D>& getCor(void) const;

  //! \brief Class constructor.
  Delaunay(void);

  /*!
    \brief Copy constructor
    \param other The Triangulation to copy
  */
  Delaunay(Delaunay const& other);

  //! \brief Default destructor.
  ~Delaunay(void);

  /*! \brief Builds the Delaunay tessellation.
    \param vp A refrence to a vector of points to be added.
    \param cpoints The edges of the processor cell.
  */
  void build_delaunay(vector<Vector2D>const& vp,vector<Vector2D> const& cpoints);

  //! \brief Dumps the Delaunay tessellation into a binary file.
  void output(void);

  /*! \brief Returns a facet.
    \param index Facet index
    \returns A reference to the selected facet.
  */
  const facet& get_facet(int index) const;

  /*! \brief Returns a coordinate of a vertice.
    \param Facet The index of the facet to check.
    \param vertice The index of the vertice in the facet.
    \param dim If dim=0 returns the x-coordinate else returns the y-coordinate.
    \returns The chosen coordinate.
  */
  double get_facet_coordinate(int Facet,int vertice, int dim);

  /*! \brief Returns a point.
    \param index The index of the point.
    \returns The chosen point.
  */
  Vector2D get_point(size_t index) const;

  /*! \brief Returns a coordinate.
    \param index The index of the point.
    \param dim If dim=0 returns the x-coordinate else returns the y-coordinate.
    \returns The chosen coordinate.
  */
  double get_cor(int index,int dim) const;

  /*! \brief Returns the number of facets.
    \returns The number of facets.
  */
  int get_num_facet(void)const;

  /*! \brief Returns the number of points
    \returns The number of points.
  */
  int get_length(void) const;

  /*! \brief Returns the last location, a number used to identify the fact that the neighbor of a facet is empty.
    \returns The last location.
  */
  int get_last_loc(void) const;

  /*! \brief Change Mesh point.
    \param index The index of the point to change.
    \param p The new point to set.
  */
  void set_point(int index, Vector2D p);

  /*! \brief Returns the area of the triangle. Negative result means the triangle isn't right handed.
    \param index The index to the facet
    \return The area
  */
  double triangle_area(int index);

  /*!
    \brief Updates the triangulation
    \param points The new set of points
	\param cpoints The points of the processor cell
  */
  void update(const vector<Vector2D>& points,vector<Vector2D> const& cpoints);

  /*!
    \brief Returns the original index of the duplicated point in Periodic Boundary conditions
    \param NewPoint The index of the duplicated point
    \return The original index
  */
  int GetOriginalIndex(int NewPoint) const;

  /*!
    \brief Returns the original length of the points (without duplicated points)
    \return The original length
  */
  int GetOriginalLength(void) const;

  /*!
    \brief Returns a refrence to the points
    \return Refrence to the points
  */
  vector<Vector2D>& GetMeshPoints(void);

  /*! \brief Returns the length of all the points (included duplicated)
    \return The length of all of the points
  */
  int GetTotalLength(void);

  /*! \brief return the facet's radius
    \param facet The facet to check
    \return The facet's radius
  */
  double GetFacetRadius(int facet) const;

  /*!
  \brief Adds a point to the cor vector. Used in periodic boundaries with AMR.
  \param vec The point to add.
  */
  void AddAditionalPoint(Vector2D const& vec);
   /*!
  \brief Gets the size of the cor vector.
  \return The size of the cor vector.
  */
  int GetCorSize(void)const;
  
  /*!
  \brief Returns the center of the circumscribed circle of a facet
  \param index The index of the facet
  \returns The circumscribed circle's center
  */
  Vector2D GetCircleCenter(int index)const;

  //! \brief Diagnostics
  delaunay_loggers::DelaunayLogger* logger;
  /*!
  \brief Builds the boundary points
  \param obc The geometrical boundary conditions
  \param edges The edges of the domain
  \return The indeces of the boundary points for each edge, can be larger than the number of edges since it include corners at the end
  */
  vector<vector<int> > BuildBoundary(OuterBoundary const* obc,vector<Edge> const& edges);
  /*!
  \brief Builds the boundary points for parallel runs
  \param obc The geometrical boundary conditions
  \param tproc The tessellation of the processors
  \param Nghost The indeces of the ghost cells (order by cpu) in the cor vector. Given as output.
  \return The indeces of the boundary points sent to each cpu and the list of cpus to talk with.
  */
  pair<vector<vector<int> >,vector<int> > BuildBoundary(OuterBoundary const* obc,Tessellation const& tproc,
	  vector<vector<int> > &Nghost);
  /*!
  \brief Adds the points to the tessellation, used for boundary points
  \param points The points to add
  */
  void AddBoundaryPoints(vector<Vector2D> const& points);
};
#endif //DELAUNAYMPI_HPP
