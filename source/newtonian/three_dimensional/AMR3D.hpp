/*! \file AMR3D.hpp
\brief Abstract class for amr in 3D
\author Elad Steinberg
*/

#ifndef AMR3D_HPP
#define AMR3D_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../common/equation_of_state.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include <boost/scoped_ptr.hpp>
#include "LinearGauss3D.hpp"
#include "hdsim_3d.hpp"

//! \brief Abstract class for cell update scheme in amr
class AMRCellUpdater3D
{
public:

	/*! \brief Calculates the computational cell
	\param extensive Extensive conserved variables
	\param eos Equation of state
	\param volume Cell volume
	\param old_cell Old computational cell
	\param tracerstickernames The names of the tracers and stickers
	\return Computational cell
	*/
	virtual ComputationalCell3D ConvertExtensiveToPrimitve3D(const Conserved3D& extensive, const EquationOfState& eos,
		double volume, ComputationalCell3D const& old_cell, TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Class destructor
	virtual ~AMRCellUpdater3D(void);
};

//! \brief Abstract class for extensive update scheme in amr
class AMRExtensiveUpdater3D
{
public:

	/*! \brief Calculates the computational cell
	\param cell Computational cell
	\param eos Equation of state
	\param volume Cell volume
	\param tracerstickernames The names of the tracers and stickers
	\return Extensive
	*/
	virtual Conserved3D ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
		double volume, TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Class destructor
	virtual ~AMRExtensiveUpdater3D(void);
};

//! \brief Simple class for extensive update scheme in amr
class SimpleAMRExtensiveUpdater3D : public AMRExtensiveUpdater3D
{
public:
	SimpleAMRExtensiveUpdater3D(void);

	Conserved3D ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
		double volume, TracerStickerNames const& tracerstickernames) const;
};

//! \brief Simple class for cell update scheme in amr
class SimpleAMRCellUpdater3D : public AMRCellUpdater3D
{
private:
	const vector<string>  toskip_;
public:
	/*!
	\brief class constructor
	\param toskip A list of sticker names to skip their cell update
	*/
	SimpleAMRCellUpdater3D(vector<string> toskip = vector<string>());

	ComputationalCell3D ConvertExtensiveToPrimitve3D(const Conserved3D& extensive, const EquationOfState& eos,
		double volume, ComputationalCell3D const& old_cell, TracerStickerNames const& tracerstickernames) const;
};

//! \brief Chooses which cells should be remove
class CellsToRemove3D
{
public:
	/*!
	\brief Finds the cells to remove
	\param tess The tesselation
	\param cells The computational cells
	\param time The sim time
	\param tracerstickernames The names of the tracers and stickers
	\return The indeces of cells to remove with a corresponding merit which decides if there are neighboring cells which one to choose to remove
	*/
	virtual std::pair<vector<size_t>, vector<double> > ToRemove(Tessellation3D const& tess,
		vector<ComputationalCell3D> const& cells, double time,
		TracerStickerNames const& tracerstickernames)const = 0;

	//! \brief Virtual destructor
	virtual ~CellsToRemove3D(void);
};

//! \brief Chooses which cells should be refined
class CellsToRefine3D
{
public:
	/*!
	\brief Finds the cells to refine
	\param tess The tesselation
	\param cells The computational cells
	\param time The sim time
	\param tracerstickernames The names of the tracers and stickers
	\return The indeces of cells to remove and the direction to split (can be given empty)
	*/
	virtual std::pair<vector<size_t>,vector<Vector3D> > ToRefine(Tessellation3D const& tess,
		vector<ComputationalCell3D> const& cells, double time, TracerStickerNames const& tracerstickernames)const = 0;

	//! \brief Virtual destructor
	virtual ~CellsToRefine3D(void);
};

//! \brief Base class for amr
class AMR3D
{
private:
	EquationOfState const& eos_;
	CellsToRefine3D const& refine_;
	CellsToRemove3D  const& remove_;
	SimpleAMRCellUpdater3D scu_;
	SimpleAMRExtensiveUpdater3D seu_;
	AMRCellUpdater3D* cu_;
	AMRExtensiveUpdater3D* eu_;
	LinearGauss3D *interp_;
	AMR3D(AMR3D const& amr);
	AMR3D& operator=(AMR3D const&);
	
#ifdef RICH_MPI

	/*! \brief Removes points because they are near the edge of a cpu domain
	\param ToRemove Candidates for AMR
	\param merits The merits for points to be removed. given as input and output. should be empty vector for refinement
	\param tess Tessellation
	\return The new indices and merits of points
	*/
	vector<size_t> RemoveNearBoundaryPoints(vector<size_t> const&ToRemove,
		Tessellation3D const& tess, vector<double> &merits)const;
#endif
	void UpdateCellsRefine(Tessellation3D &tess, vector<ComputationalCell3D> &cells, EquationOfState const& eos,
		vector<Conserved3D> &extensives, double time,
#ifdef RICH_MPI
		Tessellation3D const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;

	void UpdateCellsRemove2(Tessellation3D &tess, vector<ComputationalCell3D> &cells, vector<Conserved3D> &extensives,
		EquationOfState const& eos, double time, TracerStickerNames const& tracerstickernames
#ifdef RICH_MPI
		,Tessellation3D const& proctess
#endif
		)const;
public:
	/*!
	\brief Runs the AMR
	\param sim The sim object
	*/
	void operator() (HDSim3D &sim);


	/*! \brief Class constructor
	\param refine Refinement scheme
	\param remove Removal scheme
	\param cu Cell updater
	\param eu Extensive updater
	\param slopes Slopes
	\param eos Equation of state
	*/
	AMR3D(EquationOfState const& eos, CellsToRefine3D const& refine, CellsToRemove3D const& remove, LinearGauss3D *slopes = 0, AMRCellUpdater3D* cu = 0,
		AMRExtensiveUpdater3D* eu = 0);
	//! Class destructor
	~AMR3D();
};

#endif // AMR3D_HPP
