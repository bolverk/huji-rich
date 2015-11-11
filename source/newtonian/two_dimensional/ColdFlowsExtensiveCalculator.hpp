/*! \file ColdFlowsExtensiveCalculator.hpp
\author Elad Steinberg
\brief Cold flows extensive variable updater, uses entropy for thermal energy
*/

#ifndef COLD_EXTENSIVE_UPDATER_HPP
#define COLD_EXTENSIVE_UPDATER_HPP 1

#include "extensive_updater.hpp"
#include "simple_flux_calculator.hpp"
#include "../../misc/utils.hpp"
#include "interpolations/LinearGaussImproved.hpp"

//! \brief Cold flows extensive variable update
class ColdFlowsExtensiveCalculator : public ExtensiveUpdater
{
private:
	EquationOfState const& eos_;
	LinearGaussImproved const& interp_;
public:

  /*! \brief Class constructor
    \param eos Equation of state
    \param interp Interpolation
   */
	ColdFlowsExtensiveCalculator(EquationOfState const& eos,LinearGaussImproved const& interp);

	void operator()
		(const vector<Extensive>& fluxes,
		const PhysicalGeometry& pg,
		const Tessellation& tess,
		const double dt,
		const CacheData& cd,
		const vector<ComputationalCell>& cells,
		vector<Extensive>& extensives) const;
};

#endif // COLD_EXTENSIVE_UPDATER_HPP
