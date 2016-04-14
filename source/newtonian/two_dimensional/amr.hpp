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
#include "clipper.hpp"
#include "../test_2d/main_loop_2d.hpp"
#include "../../tessellation/polygon_overlap_area.hpp"
#include <boost/scoped_ptr.hpp>
#include "interpolations/LinearGaussImproved.hpp"

//! \brief Abstract class for cell update scheme in amr
class AMRCellUpdater
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
	virtual ComputationalCell ConvertExtensiveToPrimitve(const Extensive& extensive,const EquationOfState& eos,
		double volume,ComputationalCell const& old_cell,TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Class destructor
	virtual ~AMRCellUpdater(void);
};

//! \brief Abstract class for extensive update scheme in amr
class AMRExtensiveUpdater
{
public:

	/*! \brief Calculates the computational cell
	\param cell Computational cell
	\param eos Equation of state
	\param volume Cell volume
	\param tracerstickernames The names of the tracers and stickers
	\return Extensive
	*/
	virtual Extensive ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
		double volume,TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Class destructor
	virtual ~AMRExtensiveUpdater(void);
};

//! \brief Simple class for extensive update scheme in amr
class SimpleAMRExtensiveUpdater : public AMRExtensiveUpdater
{
public:
	Extensive ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
		double volume, TracerStickerNames const& tracerstickernames) const;
};

//! \brief Simple class for cell update scheme in amr
class SimpleAMRCellUpdater : public AMRCellUpdater
{
private:
	const vector<string>  toskip_;
public:
	/*!
	\brief class constructor
	\param toskip A list of sticker names to skip their cell update
	*/
  explicit SimpleAMRCellUpdater(vector<string> toskip);

	ComputationalCell ConvertExtensiveToPrimitve(const Extensive& extensive, const EquationOfState& eos,
		double volume, ComputationalCell const& old_cell,TracerStickerNames const& tracerstickernames) const;
};

//! \brief Chooses which cells should be remove
class CellsToRemove
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
	virtual std::pair<vector<size_t>,vector<double> > ToRemove(Tessellation const& tess,
		vector<ComputationalCell> const& cells,double time,
		TracerStickerNames const& tracerstickernames)const=0;

	//! \brief Virtual destructor
	virtual ~CellsToRemove(void);
};

//! \brief Chooses which cells should be refined
class CellsToRefine
{
public:
	/*!
	\brief Finds the cells to refine
		\param tess The tesselation
		\param cells The computational cells
		\param time The sim time
		\param tracerstickernames The names of the tracers and stickers
		\return The indeces of cells to remove
     */
	virtual vector<size_t> ToRefine(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
		TracerStickerNames const& tracerstickernames)const = 0;
	
	//! \brief Virtual destructor
	virtual ~CellsToRefine(void);
};

//! \brief Base class for amr
class AMR : public Manipulate
{
protected:

  /*! \brief Calculates the positions of the new points
    \param ToRefine points to refine
    \param tess Tessellation
    \param NewPoints output
    \param Moved displacement for periodic grid
    \param obc Outer boundary conditions
    */
#ifdef RICH_MPI
  //!    \param proc_chull Cells' convex hull
#endif // RICH_MPI
  void GetNewPoints
  (vector<size_t> const& ToRefine,
   Tessellation const& tess,
   vector<std::pair<size_t, Vector2D> > &NewPoints, 
   vector<Vector2D> &Moved,
   OuterBoundary const& obc
#ifdef RICH_MPI
	  , vector<Vector2D> const& proc_chull
#endif
	  )const;
#ifdef RICH_MPI

  /*! \brief Removes points because they are near the edge of a cpu domain
    \param ToRemove Candidates for AMR
	\param merits The merits for points to be removed. given as input and output. should be empty vector for refinement
    \param tess Tessellation
    \return The new indices and merits of points 
   */
  vector<size_t> RemoveNearBoundaryPoints(vector<size_t> const&ToRemove,
	  Tessellation const& tess,vector<double> &merits)const;
#endif

  /*! \brief Calculates the positions of the new points like AREPO
  \param ToRefine points to refine
  \param tess Tessellation
  \param NewPoints output
  \param Moved displacement for periodic grid
  \param obc Outer boundary conditions
  */
  void GetNewPoints2
	  (vector<size_t> const& ToRefine,
		  Tessellation const& tess,
		  vector<std::pair<size_t, Vector2D> > &NewPoints,
		  vector<Vector2D> &Moved,
		  OuterBoundary const& obc)const;

public:
	/*!
	\brief Runs the AMR
	\param sim The sim object
	*/
	virtual void operator() (hdsim &sim) = 0;
	/*!
	\brief Runs the refine
	\param tess The tessellation
	\param cells The computational cells
	\param eos The equation of state
	\param extensives The extensive variables
	\param time The sim time
	\param obc Outer boundary conditions
	\param tracerstickernames The names of the tracers and stickers
	*/
#ifdef RICH_MPI
  //! \param proctess Tessellation of the processes (for parallel runs)
#endif // RICH_MPI
	virtual void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
		vector<Extensive> &extensives, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const = 0;
	/*!
	\brief Runs the removal
	\param tess The tessellation
	\param cells The computational cells
	\param eos The equation of state
	\param extensives The extensive variables
	\param time The sim time
	\param obc The outer boundary conditions
	\param tracerstickernames The names of the tracers and stickers
	*/
#ifdef RICH_MPI
  //!	\param proctess Tessellation of the processes (for parallel runs)
#endif // RICH_MPI
	virtual void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const = 0;
	//! \brief Virtual destructor
	virtual ~AMR(void);
};

//! \brief Conservative amr
//! \todo Make sure AMR works with all physical geometries
class ConservativeAMR : public AMR
{
private:
	CellsToRefine const& refine_;
	CellsToRemove const& remove_;
  SimpleAMRCellUpdater scu_;
  SimpleAMRExtensiveUpdater seu_;
  const bool periodic_;
	AMRCellUpdater* cu_;
	AMRExtensiveUpdater* eu_;
	LinearGaussImproved *interp_;

	ConservativeAMR(ConservativeAMR const& other);

	ConservativeAMR& operator=(ConservativeAMR const& other);

public:
	void operator() (hdsim &sim);

  /*! \brief Class constructor
    \param refine Refinement scheme
    \param remove Removal scheme
    \param cu Cell updater
    \param eu Extensive updater
    \param slopes Slopes
	\param periodic Flag whether the boundary is periodic or not
   */
	ConservativeAMR
	(CellsToRefine const& refine,
	 CellsToRemove const& remove,
		bool periodic = false,
	 LinearGaussImproved *slopes = 0,
	 AMRCellUpdater* cu=0,
	 AMRExtensiveUpdater* eu=0);

	void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells,EquationOfState const& eos,
		vector<Extensive> &extensives,double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;

	void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos,double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;
};

//! \brief Non conservative amr
class NonConservativeAMR : public AMR
{
private:
	CellsToRefine const& refine_;
	CellsToRemove const& remove_;
  SimpleAMRCellUpdater scu_;
  SimpleAMRExtensiveUpdater seu_;
  LinearGaussImproved *interp_;
	AMRExtensiveUpdater* eu_;

	NonConservativeAMR(NonConservativeAMR const& other);

	NonConservativeAMR& operator=(NonConservativeAMR const& other);

public:
	void operator() (hdsim &sim);

  /*! \brief Class constructor
    \param refine Refinement scheme
    \param remove Removal scheme
    \param eu Extensive updater
	\param slopes Slopes
   */
  NonConservativeAMR
  (CellsToRefine const& refine,
   CellsToRemove const& remove,
	  LinearGaussImproved *slopes = 0,
   AMRExtensiveUpdater* eu = 0);

	void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
		vector<Extensive> &extensives, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;

	void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;
};

//! \brief Conservative amr using old method to split cells
//! \todo Make sure AMR works with all physical geometries
class ConservativeAMROld : public AMR
{
private:
	CellsToRefine const& refine_;
	CellsToRemove const& remove_;
	SimpleAMRCellUpdater scu_;
	SimpleAMRExtensiveUpdater seu_;
	AMRCellUpdater* cu_;
	AMRExtensiveUpdater* eu_;
	LinearGaussImproved *interp_;

	ConservativeAMROld(ConservativeAMROld const& other);

	ConservativeAMROld& operator=(ConservativeAMROld const& other);

public:
	void operator() (hdsim &sim);

	/*! \brief Class constructor
	\param refine Refinement scheme
	\param remove Removal scheme
	\param cu Cell updater
	\param eu Extensive updater
	\param slopes Slopes
	*/
	ConservativeAMROld
		(CellsToRefine const& refine,
			CellsToRemove const& remove,
			LinearGaussImproved *slopes = 0,
			AMRCellUpdater* cu = 0,
			AMRExtensiveUpdater* eu = 0);

	void UpdateCellsRefine(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
		vector<Extensive> &extensives, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;

	void UpdateCellsRemove(Tessellation &tess,
		OuterBoundary const& obc, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
		EquationOfState const& eos, double time,
#ifdef RICH_MPI
		Tessellation const& proctess,
#endif
		TracerStickerNames const& tracerstickernames)const;
};


#endif // AMR_HPP
