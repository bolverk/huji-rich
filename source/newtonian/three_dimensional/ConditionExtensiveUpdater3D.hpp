#ifndef CONDITION_EXTENSIVE_UPDATER3D_HPP
#define CONDITION_EXTENSIVE_UPDATER3D_HPP 1

#include "extensive_updater3d.hpp"
#include "LinearGauss3D.hpp"
#include "../common/equation_of_state.hpp"

//! \brief Updates the extensives based on a series of conditions and actions. Does a normal update of all cells before going into the conditions check
class ConditionExtensiveUpdater3D : public ExtensiveUpdater3D
{
public:

	//! \brief Determines the kind of cell
	class Condition3D
	{
	public:

		/*! \brief Checks if an interface satisfies a condition
		\param index The index of the cell
		\param tess Tessellation
		\param cells Computational cells
		\param time The sim time
		\param tracerstickernames The names of the tracers and stickers
		\return Whether the cell satisfies a condition.
		*/
		virtual bool operator()(size_t index,const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
			double time,TracerStickerNames const& tracerstickernames) const = 0;

		virtual ~Condition3D(void);
	};

	//! \brief Action taken to update extensive
	class Action3D
	{
	public:
		/*!
		\param fluxes Fluxes
		\param tess Tessellation
		\param dt Time step
		\param cells Computational cells
		\param extensives Extensive variables, input is after the addition of hydro fluxes
		\param index The index of the cell
		\param time The time
		\param tracerstickernames The names of the tracers and stickers
		*/
		virtual void operator()	(const vector<Conserved3D>& fluxes,const Tessellation3D& tess,const double dt,
			const vector<ComputationalCell3D>& cells,vector<Conserved3D> &extensives,size_t index,double time,
			TracerStickerNames const& tracerstickernames) const = 0;

		virtual ~Action3D(void);
	};

	/*! \brief Class constructor
	\param sequence Series of condition and action action pairs. Both have to be dynamically allocated pointers, and will be explicitly destructed upon descruction of the class
	*/
	explicit ConditionExtensiveUpdater3D(const vector<pair<const Condition3D*, const Action3D*> >& sequence);

	~ConditionExtensiveUpdater3D(void);

	void operator()(const vector<Conserved3D>& fluxes,const Tessellation3D& tess,const double dt,
		const vector<ComputationalCell3D>& cells,vector<Conserved3D>& extensives,double time,
		TracerStickerNames const& tracerstickernames) const;

private:
	const vector<pair<const Condition3D*, const Action3D*> > sequence_;
};


//! \brief Updates the extensive with entropy if needed for pressure
class ColdFlowsUpdate3D : public ConditionExtensiveUpdater3D::Action3D
{
public:

	/*! \brief Class constructor
	\param eos The equation of state
	\param ghost The ghost point generator
	\param interp The interpolation
	*/
	ColdFlowsUpdate3D(EquationOfState const& eos, Ghost3D const& ghost, LinearGauss3D const& interp);

	void operator()	(const vector<Conserved3D>& fluxes,const Tessellation3D& tess,const double dt,
		const vector<ComputationalCell3D>& cells, vector<Conserved3D> &extensives,size_t index,double time,
		TracerStickerNames const& tracerstickernames)const;
private:
	EquationOfState const& eos_;
	Ghost3D const& ghost_;
	LinearGauss3D const& interp_;
	mutable double lasttime_, dt_;
	mutable size_t entropy_index_;
	mutable boost::container::flat_map<size_t, ComputationalCell3D> ghost_cells_;
};

class ChooseAll : public ConditionExtensiveUpdater3D::Condition3D
{
public:

	bool operator()(size_t /*index*/, const Tessellation3D& /*tess*/, const vector<ComputationalCell3D>& /*cells*/,
		double /*time*/, TracerStickerNames const& /*tracerstickernames*/) const
	{
		return true;
	}
};

class RegularExtensiveUpdate3D : public ConditionExtensiveUpdater3D::Action3D
{
public:


	void operator()	(const vector<Conserved3D>& fluxes, const Tessellation3D& tess, const double dt,
		const vector<ComputationalCell3D>& cells, vector<Conserved3D> &extensives, size_t index, double time,
		TracerStickerNames const& tracerstickernames)const;
};

#endif // CONDITION_EXTENSIVE_UPDATER3D_HPP
