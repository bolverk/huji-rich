/*! \file RemovalStrategy.hpp
  \brief Abstract class for coarsening scheme
  \author Elad Steinberg
 */

#ifndef REMOVALSTRATEGY_HPP
#define REMOVALSTRATEGY_HPP 1

#include "../../tessellation/tessellation.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"
#include <algorithm>
#include "../../misc/universal_error.hpp"

/*! \brief Abstract class for derefinment strategies
	\author Elad Steinberg
*/
class RemovalStrategy
{
public:
	/*!
		\brief Removal abstract class. Can't remove neighboring cells or cells near peiodic boundary. Use CheckOutput to check correctness
		\param tess The tessellation
		\param cells The hydro primitives
		\param time The sim time
		\return The cells to remove
	*/
	virtual vector<int> CellsToRemove
	(Tessellation const& tess,
	 vector<ComputationalCell> const& cells,
	 double time)const=0;
	
	/*!
		\brief Checks if the removed list is good, throws an error if not
		\param tess The tessellation
		\param ToRemove The list of points to remove
	*/
	void CheckOutput(Tessellation const& tess,vector<int> & ToRemove) const;

	/*! \brief Removes neighboring points and cells near the boundary which are not rigid walls
	\param merits The vector of merits that decides which one of the neighbors to keep (the one with the higher merit)
	\param candidates The list of points to remove, assumed to be sorted
	\param tess The tessellation
	\return The list of points to remove without neighboring points
	*/
	vector<int> RemoveNeighbors(vector<double> const& merits,vector<int> const&
		candidates,Tessellation const& tess) const;
	/*! \brief Removed from the list cells near periodic boundaries
	\param ToRemove List of candidates
	\param tess The tessellation
	\return The corrected list
	*/
	vector<int> RemoveNearBoundary(vector<int> const& ToRemove,Tessellation
		const& tess)const;

	//! \brief Virtual destructor
  virtual ~RemovalStrategy(void);
};

#endif //REMOVALSTRATEGY_HPP
