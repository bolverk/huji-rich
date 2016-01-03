/*!
\file AreaFix.hpp
\brief Fixes the area inconsistency problem
\author Elad Steinberg
*/

#ifndef AREAFIX_HPP
#define AREAFIX_HPP 1

#include "../../tessellation/tessellation.hpp"
#include "../common/equation_of_state.hpp"
#include "computational_cell_2d.hpp"
#include "extensive.hpp"
#include "OuterBoundary.hpp"
#include "../../misc/utils.hpp"
#include "../../tessellation/geotests.hpp"

/*!
\brief Fixes the flux to propely converge second order
\param tessold The tessellation at the beginning of the time step
\param tessmid The tessellation at the mid time step (old for first fix)
\param tessnew The tessellation at the end of the time step (mid for first fix)
\param pointvelocity The mesh point velocities
\param dt The time step
\param cells The computational cells
\param fluxes The fluxes
\param fv The edge velocities
\param outer The outer boundary conditions
\param eos The equation of state
\return The fix for each cells
*/
vector<Extensive> FluxFix2(Tessellation const& tessold, Tessellation const& tessmid,
	Tessellation const& tessnew, vector<Vector2D> const& pointvelocity, double dt,
	vector<ComputationalCell> const& cells, vector<Extensive> const& fluxes,
	vector<Vector2D> const& fv, OuterBoundary const& outer, EquationOfState const& eos);

#endif //AREAFIX_HPP
