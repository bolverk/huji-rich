/*! \file amr.hpp
\brief Abstract class for amr
\author Elad Steinberg
*/

#ifndef AMR_HPP
#define AMR_HPP 1

#include "computational_cell_2d.hpp"
#include "extensive.hpp"
#include "../common/equation_of_state.hpp"
#include "OuterBoundary.hpp"
#include "../../tessellation/tessellation.hpp"
#include "../../tessellation/ConvexHull.hpp"
#include "../../clipper/clipper.hpp"
#include "../test_2d/main_loop_2d.hpp"
#include "../../tessellation/polygon_overlap_area.hpp"
#include <boost/scoped_ptr.hpp>

//! \brief Abstract class for cell update scheme in amr
class AMRCellUpdater
{
public:

	/*! \brief Calculates the computational cell
	\param intensive Intensive conserved variables (per volume)
	\param eos Equation of state
	\return Computational cell
	*/
	virtual ComputationalCell ConvertExtensiveToPrimitve(const Extensive& extensive,const EquationOfState& eos,
		double volume,ComputationalCell const& old_cell) const = 0;

	//! \brief Class destructor
	virtual ~AMRCellUpdater(void);
};

// !\brief Abstract class for extensive update scheme in amr
class AMRExtensiveUpdater
{
public:

	/*! \brief Calculates the computational cell
	\param intensive Intensive conserved variables (per volume)
	\param eos Equation of state
	\return Computational cell
	*/
	virtual Extensive ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
		double volume) const = 0;

	//! \brief Class destructor
	virtual ~AMRExtensiveUpdater(void);
};

// !\brief Simple class for extensive update scheme in amr
class SimpleAMRExtensiveUpdater : public AMRExtensiveUpdater
{
public:
	Extensive ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
		double volume) const;
};

// !\brief Simple class for cell update scheme in amr
class SimpleAMRCellUpdater : public AMRCellUpdater
{
public:
	ComputationalCell ConvertExtensiveToPrimitve(const Extensive& extensive, const EquationOfState& eos,
		double volume, ComputationalCell const& old_cell) const;
};


class CellsToRemove
{
public:
	virtual std::pair<vector<size_t>,vector<double> > ToRemove(Tessellation const& tess,
		vector<ComputationalCell> const& cells,double time)const=0;

	virtual ~CellsToRemove(void);
};

class CellsToRefine
{
public:
	virtual vector<size_t> ToRefine(Tessellation const& tess, vector<ComputationalCell> const& cells, double time)const = 0;
	
	virtual ~CellsToRefine(void);
};


class AMR : public Manipulate
{
public:
	virtual void operator() (hdsim &sim) = 0;

	virtual void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
		vector<Extensive> &extensives, double time)const = 0;

	virtual void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos, double time)const = 0;

	void GetNewPoints(vector<size_t> const& ToRefine, Tessellation const& tess,
		vector<std::pair<size_t, Vector2D> > &NewPoints, vector<Vector2D> &Moved,
		OuterBoundary const& obc)const;

	virtual ~AMR(void);
};

//! \todo Make sure AMR works with all physical geometries
class ConservativeAMR : public AMR
{
private:
	CellsToRefine const& refine_;
	CellsToRemove const& remove_;
	static SimpleAMRCellUpdater scu_;
	static SimpleAMRExtensiveUpdater seu_;
	AMRCellUpdater const& cu_;
	AMRExtensiveUpdater const& eu_;

	vector<size_t> RemoveNeighbors(vector<double> const& merits, vector<size_t> const&
		candidates, Tessellation const& tess) const;

public:
	void operator() (hdsim &sim);

	ConservativeAMR(CellsToRefine const& refine, CellsToRemove const& remove, AMRCellUpdater const& cu = scu_,
		AMRExtensiveUpdater const& eu = seu_);

	void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells,EquationOfState const& eos,
		vector<Extensive> &extensives,double time)const;

	void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos,double time)const;
};

class NonConservativeAMR : public AMR
{
private:
	CellsToRefine const& refine_;
	CellsToRemove const& remove_;
	static SimpleAMRCellUpdater scu_;
	static SimpleAMRExtensiveUpdater seu_;
	AMRCellUpdater const& cu_;
	AMRExtensiveUpdater const& eu_;

public:
	void operator() (hdsim &sim);

	NonConservativeAMR(CellsToRefine const& refine, CellsToRemove const& remove, AMRCellUpdater const& cu = scu_,
		AMRExtensiveUpdater const& eu = seu_);

	void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
		vector<Extensive> &extensives, double time)const;

	void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos, double time)const;
};

#endif // AMR_HPP
