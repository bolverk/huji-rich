/*! \file GhostPointGenerator.hpp
\brief Abstract class for creating computationalcells of ghost points
\author Elad Steinberg
*/

#ifndef GHOST_POINT_GENERATOR_HPP
#define GHOST_POINT_GENERATOR_HPP 1

#include "../common/hydrodynamic_variables.hpp"
#include "../../tessellation/tessellation.hpp"
#include "computational_cell_2d.hpp"
#include <map>

/*! \brief Abstract class for creating ghost points
\author Elad Steinberg
*/
class GhostPointGenerator
{
public:
	/*!
	\brief Calculates the ghost points
	\param tess The tessellation
	\param cells The computational cells
	\return A map where the key is the index of the ghost cell and the value is its' comuptational cell
	*/
	virtual std::map<size_t,ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells) const = 0;
	//! \brief Virtual destructor
	virtual ~GhostPointGenerator(void);
	/*!
	\brief Finds the indeces of the outer edges points
	\param tess The tessellation
	\return The indeces of the outer edges and whether the ghost is the first neighbor (1) or the second (2)
	*/
	vector<std::pair<size_t,size_t> > GetOuterEdgesIndeces(Tessellation const& tess)const;
};

#endif // GHOST_POINT_GENERATOR_HPP
