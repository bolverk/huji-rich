/*! \file RefineStrategy.hpp
  \brief Abstract class for refinement scheme
  \author Elad Steinberg
 */

#ifndef REFINESTRATEGY_HPP
#define REFINESTRATEGY_HPP 1

#include <algorithm>
#include "../../tessellation/tessellation.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"
#include "../../misc/utils.hpp"
/*! \brief Abstract class for refinment strategies
	\author Elad Steinberg
*/
class RefineStrategy
{
public:
	//! \brief Class constructor
  RefineStrategy(void);
  /*!
  \brief Calculates the cells to be refined
  \param tess The tessellation
  \param cells The primitive cells
  \param time The simulation time
  \param directions The directions to move the splitted points, can be given empty
  \param Removed A list of the cells that were removed in the last cell removal
  \return A list of the cells to refine
  */
	virtual vector<int> CellsToRefine
	(Tessellation const& tess,
	 vector<ComputationalCell> const& cells,
	 double time,
	 vector<Vector2D> &directions,
	 vector<int> const& Removed)=0;
	/*!
	\brief Removes cells that were splitted in the last time step
	\param ToRefine The list of candidate cells to split
	\param Npoints The number of points in the tessellation
	\param directions The directions to move the splitted points, can be given empty
	\param Removed A list of the cells that were removed in the last cell removal
	\param tess The tessellation
	\return A new list, where the points that were split last time step are removed from
	*/
	vector<int> RemoveDuplicatedLately(vector<int> const& ToRefine,
		int Npoints,vector<Vector2D> &directions,vector<int> const& Removed,
		Tessellation const& tess);

	/*! \brief Removed from the list cells near periodic boundaries
	\param ToRefine List of candidates
	\param directions The split directions, can be given empty
	\param tess The tessellation
	\return The corrected list
	*/
	vector<int> RemoveNearBoundary(vector<int> const& ToRefine,vector<Vector2D>
		&directions,Tessellation const& tess);
	//! \brief The cells that were refined in the previous time step
	vector<int> refined_old;
//! \brief Virtual destructor
  virtual ~RefineStrategy(void);
};

/*! \brief Finds the best way to spit the cell by finding the line connecting the mesh point and the cell's CM
\param tess The tessellation
\param PointToRefine The index of the refined cell
\param edges The edges of the refined cell
\param R The width of the refined cell
\param normal The normal to the split direction, this is one of the outputs
\return slope The direction where the new point should be at
*/
Vector2D FindBestSplit(Tessellation const* tess,int PointToRefine,
	vector<Edge> const& edges,double R,Vector2D &normal);
#endif // REFINESTRATEGY_HPP
