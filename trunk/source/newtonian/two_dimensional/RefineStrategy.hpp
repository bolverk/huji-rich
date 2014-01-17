/*! \file RefineStrategy.hpp
  \brief Abstract class for refinement scheme
  \author Elad Steinberg
 */

#ifndef REFINESTRATEGY_HPP
#define REFINESTRATEGY_HPP 1

#include <algorithm>
#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
#include "../../misc/utils.hpp"
#include <iostream>
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
  \param tracers The tracers
  \param time The simulation time
  \param directions The directions to move the splitted points, can be given empty
  \param Removed A list of the cells that were removed in the last cell removal
  \return A list of the cells to refine
  */
	virtual vector<int> CellsToRefine(Tessellation const* tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double time,vector<Vector2D> &directions,vector<int> const& Removed)=0;
	/*!
	\brief Removes cells that were splitted in the alst time step
	\param ToRefine The list of candidate cells to split
	\param Npoints The number of points in the tessellation
	\param directions The directions to move the splitted points, can be given empty
	\param Removed A list of the cells that were removed in the last cell removal
	\return A new list, where the points that were split last time step are removed from
	*/
	vector<int> RemoveDuplicatedLately(vector<int> const& ToRefine,
		int Npoints,vector<Vector2D> &directions,vector<int> const& Removed);
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
