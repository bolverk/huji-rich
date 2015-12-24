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

vector<Extensive> FluxFix2(Tessellation const& tessold, Tessellation const& tessmid,
	Tessellation const& tessnew, vector<Vector2D> const& pointvelocity, double dt,
	vector<ComputationalCell> const& cells, vector<Extensive> const& fluxes,
	vector<Vector2D> const& fv, OuterBoundary const& outer, EquationOfState const& eos);

#endif //AREAFIX_HPP
