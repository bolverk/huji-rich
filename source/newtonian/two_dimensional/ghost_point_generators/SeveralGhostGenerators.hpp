/*! \file SeveralGhostGenerators.hpp
\brief Class for creating computationalcells of ghost points from several different methods
\author Elad Steinberg
*/

#ifndef SEVERAL_GHOST_GENERATOR_HPP
#define SEVERAL_GHOST_GENERATOR_HPP 1

#include "../GhostPointGenerator.hpp"

/*! \brief Class for chhosing which ghost generator to use
\author Elad Steinberg
*/
class GhostCriteria
{
public:
	/*!
	\brief Chooses the ghost generator
	\param tess The tessellation
	\param index The index of the mesh point to calculate for
	\return The index in the ghost generator vector to choose from
	*/
	virtual size_t GhostChoose(Tessellation const& tess, int index)const = 0;

  virtual ~GhostCriteria(void);
};

/*! \brief Class for creating computationalcells of ghost points from several different methods
\author Elad Steinberg
*/
class SeveralGhostGenerators : public GhostPointGenerator
{
private:
	vector<GhostPointGenerator*> ghosts_;
	GhostCriteria const& ghost_chooser_;
public:

  /*! \brief Class constructor
    \param ghosts List of ghost generators
    \param ghostchooser Criteria for when to use each ghost generator
   */
	SeveralGhostGenerators(vector<GhostPointGenerator*> ghosts,GhostCriteria const& ghostchooser);

	boost::container::flat_map<size_t, ComputationalCell> operator() (const Tessellation& tess,
		const vector<ComputationalCell>& cells,double time,TracerStickerNames const&
		tracerstickernames) const;

	Slope GetGhostGradient(const Tessellation& tess,
		const vector<ComputationalCell>& cells, const vector<Slope>& gradients,
		size_t ghost_index, double time, const Edge& edge, TracerStickerNames const&
		tracerstickernames) const;
};

#endif // SEVERAL_GHOST_GENERATOR_HPP
