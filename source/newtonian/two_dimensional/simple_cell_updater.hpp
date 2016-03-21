/*! \file simple_cell_updater.hpp
  \author Almog Yalinewich
  \brief Simple cell updater
*/

#ifndef SIMPLE_CELL_UPDATER_HPP
#define SIMPLE_CELL_UPDATER_HPP 1

#include "cell_updater_2d.hpp"
#include "../../misc/utils.hpp"

using std::vector;

//! \brief Simple cell updater
class SimpleCellUpdater : public CellUpdater
{
public:

	//! \brief Abstract class to determine cell type
	class Condition
	{
	public:

		virtual ~Condition(void);

		/*! \brief Checks if a cell satisfies a certain condition
		  \param tess Tessellation
		  \param pg Physical geometry
		  \param eos Equation of state
		  \param extensives Extensives
		  \param cells Computational cells
		  \param cd Cached data
		  \param index Cell index
		  \param tracerstickernames THe names of the stickers and tracers
		  \return True if condition is met, false otherwise
		 */
		virtual bool operator()
			(const Tessellation& tess,
				const PhysicalGeometry& pg,
				const EquationOfState& eos,
				const vector<Extensive>& extensives,
				const vector<ComputationalCell>& cells,
				const CacheData& cd,
				const size_t index,
				TracerStickerNames const& tracerstickernames) const = 0;
	};

	//! \brief Action taken to calculate cell
	class Action
	{
	public:

		virtual ~Action(void);

		/*! \brief Calculates the value of the primitive variables in a computational cell
		  \param tess Tessellation
		  \param pg Physical geometry
		  \param eos Equation of state
		  \param extensives Extensive variables
		  \param cells Computational cells
		  \param cd Cached data
		  \param index Cell index
		  \param tracerstickernames THe names of the stickers and tracers
		  \return Computational cell
		 */
		virtual ComputationalCell operator()
			(const Tessellation& tess,
				const PhysicalGeometry& pg,
				const EquationOfState& eos,
				const vector<Extensive>& extensives,
				const vector<ComputationalCell>& cells,
				const CacheData& cd,
				const size_t index,
				TracerStickerNames const& tracerstickernames)const = 0;
	};

	/*! \brief Class constructor
	  \param sequence List of rules for cells that are calculated in a special way
	 */
	explicit SimpleCellUpdater
		(const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence =
			vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> >());

	vector<ComputationalCell> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			vector<Extensive>& extensives,
			const vector<ComputationalCell>& old,
			const CacheData& cd,
			TracerStickerNames const& tracerstickernames) const;

	~SimpleCellUpdater(void);

private:
	const vector<pair<const SimpleCellUpdater::Condition*, const SimpleCellUpdater::Action*> > sequence_;
	const string entropy_;
};

//! \brief Checks if a cell contains a certain sticker
class HasSticker : public SimpleCellUpdater::Condition
{
public:

	/*! \brief Class constructor
	  \param sticker_name Sticker name
	 */
	explicit HasSticker(const string& sticker_name);

	bool operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			const vector<Extensive>& extensives,
			const vector<ComputationalCell>& cells,
			const CacheData& cd,
			const size_t index, 
			TracerStickerNames const& tracerstickernames) const;

private:
	const string sticker_name_;
};

//! \brief Prevents certain cells from being updated
class SkipUpdate : public SimpleCellUpdater::Action
{
public:

	SkipUpdate(void);

	ComputationalCell operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const EquationOfState& eos,
			const vector<Extensive>& extensives,
			const vector<ComputationalCell>& cells,
			const CacheData& cd,
			const size_t index,
			TracerStickerNames const& tracerstickernames) const;
};

#endif // SIMPLE_CELL_UPDATER_HPP
