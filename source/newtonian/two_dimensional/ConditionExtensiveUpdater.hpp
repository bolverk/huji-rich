#ifndef CONDITION_EXTENSIVE_UPDATER_HPP
#define CONDITION_EXTENSIVE_UPDATER_HPP 1

#include "extensive_updater.hpp"
#include "interpolations/LinearGaussImproved.hpp"

//! \brief Updates the extensives based on a series of conditions and actions
class ConditionExtensiveUpdater : public ExtensiveUpdater
{
public:

	//! \brief Determines the kind of cell
	class Condition
	{
	public:

		/*! \brief Checks if an interface satisfies a condition
		\param index The index of the cell
		\param tess Tessellation
		\param cells Computational cells
		\param time The sim time
		\return Whether the cell satisfies a condition.
		*/
		virtual bool operator()
			(size_t index,
				const Tessellation& tess,
				const vector<ComputationalCell>& cells,
				double time,
				TracerStickerNames const& tracerstickernames) const = 0;

		virtual ~Condition(void);
	};

	//! \brief Action taken to update extensive
	class Action
	{
	public:
		/*!
		\param fluxes Fluxes
			\param pg Physical geometry
			\param tess Tessellation
			\param dt Time step
			\param cd Cache data
			\param cells Computational cells
			\param extensive Extensive variable, input is after the addition of hydro fluxes
			\param index The index of the cell
			\param time The time
		*/
		virtual void operator()
			(const vector<Extensive>& fluxes,
				const PhysicalGeometry& pg,
				const Tessellation& tess,
				const double dt,
				const CacheData& cd,
				const vector<ComputationalCell>& cells,
				Extensive& extensive,
				size_t index,
				double time,
				TracerStickerNames const& tracerstickernames) const = 0;

		virtual ~Action(void);
	};

	/*! \brief Class constructor
	\param sequence Series of condition and action action pairs. Both have to be dynamically allocated pointers, and will be explicitly destructed upon descruction of the class
	*/
	explicit ConditionExtensiveUpdater
		(const vector<pair<const Condition*, const Action*> >& sequence);

	~ConditionExtensiveUpdater(void);

	void operator()(const vector<Extensive>& fluxes,
		const PhysicalGeometry& pg,
		const Tessellation& tess,
		const double dt,
		const CacheData& cd,
		const vector<ComputationalCell>& cells,
		vector<Extensive>& extensives,
		double time,
		TracerStickerNames const& tracerstickernames) const;

private:
	const vector<pair<const Condition*, const Action*> > sequence_;
};


//! \brief Updates the extensive with entropy if needed for pressure
class ColdFlowsUpdate : public ConditionExtensiveUpdater::Action
{
public:

	/*! \brief Class constructor
	\param eos The equation of state
	\param interp The interpolation
	*/
	ColdFlowsUpdate(EquationOfState const& eos, LinearGaussImproved const& interp);

	void operator()
		(const vector<Extensive>& fluxes,
			const PhysicalGeometry& pg,
			const Tessellation& tess,
			const double dt,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			Extensive& extensive,
			size_t index,
			double time, 
			TracerStickerNames const& tracerstickernames)const;
private:
	EquationOfState const& eos_;
	LinearGaussImproved const& interp_;
};

#endif // CONDITION_EXTENSIVE_UPDATER_HPP
