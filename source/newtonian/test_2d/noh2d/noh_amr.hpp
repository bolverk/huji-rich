/*! \file noh_amr.hpp
  \brief AMR scheme for the noh problem
  \author Elad Steinberg
 */

#ifndef NOH_AMR_HPP
#define NOH_AMR_HPP 1

#include "../../../../source/newtonian/two_dimensional/RefineStrategy.hpp"
#include "../../../../source/newtonian/two_dimensional/RemovalStrategy.hpp"
//! \brief Noh problem Refinement strategy class
class NohRefine : public RefineStrategy
{
private:
	double Vmax_;
public:
	/*!
	\brief Class constructor
	\param Vmax The volume above to refine the cell
	*/
  explicit NohRefine(double Vmax);
	//! \brief Class destructor
	~NohRefine();

	vector<int> CellsToRefine
	(Tessellation const& tess,
	 vector<ComputationalCell> const& cells,
	 double time,
	 vector<Vector2D> &directions,
	 vector<int> const& Removed);
};

//! \brief Noh problem Removal strategy class
class NohRemove : public RemovalStrategy
{
private:
	double Vmin_;
public:
	/*!
	\brief Class constructor
	\param Vmin The minimum volume below which to remove the cell
	*/
  explicit NohRemove(double Vmin);
	//! \brief Class destructor
	~NohRemove();

	vector<int> CellsToRemove
	(Tessellation const& tess,
	 vector<ComputationalCell> const& cells, 
	 double time)const;
};

#endif //NOH_AMR_HPP
