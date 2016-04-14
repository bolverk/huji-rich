/*! \file tessellation.hpp
  \brief Abstract class for the tessellation
  \author Elad Steinberg
 */

#ifndef TESSELLATION_HPP
#define TESSELLATION_HPP 1

#include <vector>
#include "geometry.hpp"

#include "Edge.hpp"

using std::vector;

class OuterBoundary;

/*! \brief Abstract class for tessellation
  \author Elad Steinberg
*/
class Tessellation
{
public:
/*!
	\brief Calculates the velocity of a single edge
	\param wl The velocity of the left mesh point
	\param wr The velocity of the right mesh point
	\param rL The location of the left mesh point
	\param rR The location of the right mesh point
	\param f The centroid of the edge
	\returns The edge's velocity
*/
  virtual Vector2D CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL,
				    Vector2D rR,Vector2D f) const=0;

  /*! \brief Initialises the tessellation
    \param points Initial position of mesh generating points
    \param bc Boundary conditions of the computational domain
	\param HilbertOrder Should the points be rearranged before insertion
   */
  virtual void Initialise(vector<Vector2D> const& points, OuterBoundary const* bc, bool HilbertOrder = true) = 0;

#ifdef RICH_MPI
  /*! \brief Initialises the tessellation
  \param points Initial position of mesh generating points
  \param tess The tessellation of the processors
  \param outer The geometric outer boundary conditions
  \param HilbertOrder Should the points be rearranged before insertion
  */
  virtual void Initialise(vector<Vector2D> const& points,Tessellation const& tess,
	  OuterBoundary const* outer, bool HilbertOrder = true) = 0;
#endif

  /*!
  \brief Update the tessellation
  \param points The new positions of the mesh generating points
  \param HilbertOrder Should the points be rearranged before insertion
  \return The indeces of sort (if done, else empty)
   */
  virtual vector<int> Update(const vector<Vector2D>& points,bool HilbertOrder=false) = 0;

#ifdef RICH_MPI
  /*!
  \brief Update the tessellation
  \param points The new positions of the mesh generating points
  \param tess The tessellation of the processors
  \param HOrder Should the points be rearranged before insertion
  \return The indeces of sort (if done, else empty)
   */
  virtual vector<int> Update(const vector<Vector2D>& points, const Tessellation& tess, bool HOrder = false) = 0;
#endif // RICH_MPI

  /*! \brief Get Total number of mesh generating points
    \return Number of mesh generating points
   */
  virtual int GetPointNo(void) const = 0;

  /*! \brief Returns Position of mesh generating point
    \param index Mesh generating point index
    \return Position of mesh generating point
   */
  virtual Vector2D GetMeshPoint(int index) const = 0;

  /*! \brief Returns Position of Cell's CM
    \param index Mesh generating point index (the cell's index)
    \return Position of CM
   */
  virtual Vector2D const& GetCellCM(int index) const = 0;

  /*! \brief Returns the total number of faces
    \return Total number of faces
   */
  virtual int GetTotalSidesNumber(void) const = 0;

  /*! \brief Returns reference to the list of all edges
    \return Reference to the list of all edges
   */
  virtual const vector<Edge>& getAllEdges(void) const = 0; 

  /*! \brief Returns edge (interface between cells)
    \param index Face index
    \return Interface between cells
   */
  virtual Edge const& GetEdge(int index) const= 0;

  /*! \brief Returns the effective width of a cell
    \param index Cell index
    \return Effective cell width
   */
  virtual double GetWidth(int index) const = 0;

  /*! \brief Returns the volume of a cell
    \param index Cell index
    \return Cell volume
   */
  virtual double GetVolume(int index) const = 0;

  /*! \brief Returns the indexes of a cell's edges
    \param index Cell index
    \return Cell edges
   */
  virtual vector<int>const& GetCellEdges(int index) const =0;

/*!
	\brief Returns the original index of the duplicated point
	\param point The index of the duplicated point
	\return The original point index
*/
virtual int GetOriginalIndex(int point) const;

/*!
	\brief Returns a reference to the point vector
	\returns The reference
*/
virtual vector<Vector2D>& GetMeshPoints(void)=0;

/*!
	\brief Returns the indeces of the neighbors
	\param index The cell to check
	\return The neighbors
*/
virtual vector<int> GetNeighbors(int index)const=0;

/*!
\brief Returns the indeces of the neighbors
\param index The cell to check
\param neigh the indeces of the neighbors given as output
*/
virtual void GetNeighbors(int index,vector<int> &neigh)const = 0;

/*!
	\brief Cloning function
*/
virtual Tessellation* clone(void) const=0;

virtual ~Tessellation(void)=0;

/*! \brief Returns if the cell is adjacent to a boundary
	\param index The cell to check
	\return If near boundary
*/
virtual bool NearBoundary(int index) const = 0;

	/*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points
	\return The sent points
	*/
  virtual vector<vector<int> >& GetDuplicatedPoints(void)=0;
  /*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points
	\return The sent points
	*/
  virtual vector<vector<int> >const& GetDuplicatedPoints(void)const=0;
  	/*!
	\brief Returns the indeces of the processors with whom ghost points where exchanged
	\return The list of processors
	*/
  virtual vector<int> GetDuplicatedProcs(void)const=0;
  	/*!
	\brief Returns the indeces of the points that where sent to other processors
	\return The sent points
	*/
  virtual vector<vector<int> >const& GetSentPoints(void)const=0;
  	/*!
	\brief Returns the indeces of the processors with whom points where exchanged
	\return The list of processors
	*/
  virtual vector<int> GetSentProcs(void)const=0;

  /*!
  \brief Returns the indeces of the points that remain with the processor after the ne processor mesh is built
  \return The indeces of the points
  */
  virtual vector<size_t> GetSelfPoint(void)const=0;

  /*!
  \brief Returns the indeces of each ghost point in the vector of points that the tessellation holds
  \return The indeces where the outer index is the index of the sent processor
  */
  virtual vector<vector<int> >& GetGhostIndeces(void)=0;
  /*!
  \brief Returns the indeces of each ghost point in the vector of points that the tessellation holds
  \return The indeces where the outer index is the index of the sent processor
  */
  virtual vector<vector<int> >const& GetGhostIndeces(void)const=0;
  /*!
  \brief Returns the total number of points (including ghost)
  \return The total number of points
  */
  virtual int GetTotalPointNumber(void)const=0;

  /*!
  \brief Returns the center of masses of the cells
  \return The CM's
  */
  virtual vector<Vector2D>& GetAllCM(void)=0;

  /*! \brief Retrieves vicarious neighbors
    \param result Output
    \param point Mesh generating point index
   */
  virtual void GetNeighborNeighbors(vector<int> &result, int point)const=0;
};
#endif // TESSELLATION_HPP
