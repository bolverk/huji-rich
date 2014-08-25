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
class HydroBoundaryConditions;

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
		Vector2D rR,Vector2D f)const=0;
/*!
	\brief Calculates the velocities of the edges
	\param hbc The hydro boundary conditions
	\param point_velocities The velocities of the mesh points
	\param time The sim time
	\returns The velocities of the edges
*/
virtual vector<Vector2D> calc_edge_velocities(HydroBoundaryConditions const* hbc,
	vector<Vector2D> const& point_velocities,double time)const=0;

  /*! \brief Initialises the tessellation
    \param points Initial position of mesh generating points
    \param bc Boundary conditions of the computational domain
   */
  virtual void Initialise(vector<Vector2D> const& points,OuterBoundary const* bc) = 0;

  /*! \brief Initialises the tessellation
    \param points Initial position of mesh generating points
    \param tess The tessellation of the processors
	\param outer The geometric outer boundary conditions
   */
  virtual void Initialise(vector<Vector2D> const& points,Tessellation const& tess,
	  OuterBoundary const* outer) = 0;

  /*!
  \brief Update the tessellation
  \param points The new positions of the mesh generating points
   */
  virtual void Update(vector<Vector2D> const& points) = 0;

  /*!
  \brief Update the tessellation
  \param points The new positions of the mesh generating points
  \param tess The tessellation of the processors
   */
  virtual void Update(vector<Vector2D> const& points,Tessellation const& tess) = 0;

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
	\brief Returns the indeces of the neighbors not including -1
	\param index The cell to check
	\return The neighbors
*/

virtual vector<int> GetNeighbors(int index)const=0;

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
	\brief An efficient method of removing cells from the tessellation. Currently works only for cells not near a periodic boundary. Doesn't work with self gravity yet
	\param ToRemovevector A list of the cells to remove. There SHOULDN'T be any neighboring cells in the list
	\param VolIndex The index in the old tessellation of the cells that had their volume changed
	\param Volratio The ratio of the removed cell's volume that went into the cell from VolIndex
*/
virtual void RemoveCells(vector<int> &ToRemovevector,vector<vector<int> > &VolIndex,
	vector<vector<double> > &Volratio)=0;

/*!
	\brief Refine cells by splitting them
	\param ToRefine The indexes of the cells to refine. SHOULDN'T be near periodic boundary
	\param directions Vector between the two new points
	\param alpha How close (in cell radius) to place the new mesh point
*/
virtual void RefineCells(vector<int> const& ToRefine,vector<Vector2D> const&
	directions,double alpha)=0;

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
  virtual vector<int> GetSelfPoint(void)const=0;

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
};

#endif // TESSELLATION_HPP