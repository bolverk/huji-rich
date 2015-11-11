/*! \file ColdFlowsCellUpdate.hpp
\author Elad Steinberg
\brief Cell updater for cold flows
*/

#ifndef COLDFLOWSCELLUPDATER_HPP
#define COLDFLOWSCELLUPDATER_HPP 1

#include "simple_cell_updater.hpp"

using std::vector;

//! \brief Simple cell updater
class ColdFlowsCellUpdate : public CellUpdater
{
public:

	/*! \brief Class constructor
	\param sequence List of rules for cells that are calculated in a special way
	*/
	ColdFlowsCellUpdate
		(const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence =
			vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >());

	vector<ComputationalCell> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			vector<Extensive>& extensives,
			const vector<ComputationalCell>& old,
			const CacheData& cd) const;

private:
	SimpleCellUpdater scu_;
};


#endif // COLDFLOWSCELLUPDATER_HPP
