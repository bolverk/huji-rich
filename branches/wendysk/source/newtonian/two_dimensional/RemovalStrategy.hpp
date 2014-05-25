#ifndef REMOVALSTRATEGY_HPP
#define REMOVALSTRATEGY_HPP 1

#include "../../tessellation/tessellation.hpp"
#include "../common/hydrodynamic_variables.hpp"
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
		\param tracers The primitve tracers
		\param time The sim time
		\return The cells to remove
	*/
	virtual vector<int> CellsToRemove(Tessellation const* tess,
		vector<Primitive> const& cells,vector<vector<double> > const& tracers,
		double time)const=0;
	/*!
		\brief Checks if the removed list is good, throws an error if not
		\param tess The tessellation
		\param ToRemove The list of points to remove
	*/
	void CheckOutput(Tessellation const* tess,vector<int> & ToRemove)const;

	/*! \brief Removes neighboring points
	\param merits The vector of merits that decides which one of the neighbors to keep (the one with the higher merit)
	\param candidates The list of points to remove, assumed to be sorted
	\param tess The tessellation
	\return The list of points to remove without neighboring points
	*/
	vector<int> RemoveNeighbors(vector<double> const& merits,vector<int> const& 
		candidates,Tessellation const* tess) const;

	//! \brief Virtual destructor
  virtual ~RemovalStrategy(void);
};

#endif //REMOVALSTRATEGY_HPP
