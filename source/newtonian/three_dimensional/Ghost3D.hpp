/*! \file Ghost3D.hpp
\brief Abstract class for creating computationalcells of ghost points
\author Elad Steinberg
*/

#ifndef GHOST_POINT_GENERATOR_HPP
#define GHOST_POINT_GENERATOR_HPP 1

#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "computational_cell.hpp"
#include <boost/container/flat_map.hpp>

/*! \brief Abstract class for creating ghost points
\author Elad Steinberg
*/
class Ghost3D
{
public:
	/*!
	\brief Calculates the ghost points
	\param tess The tessellation
	\param cells The computational cells
	\param time The time
	\param tracerstickernames The names of the tracers and stickers
	\param res A map where the key is the index of the ghost cell and the value is its' comuptational cell
	*/
	virtual void operator() (const Tessellation3D& tess,
		const vector<ComputationalCell3D>& cells, double time, TracerStickerNames const&
		tracerstickernames, boost::container::flat_map<size_t, ComputationalCell3D> &res) const = 0;

	/*!
	\brief Calculates the gradients for the ghost cells
	\param tess The tessellation
	\param cells The computational cells
	\param gradients The gradients for the non-ghost cells
	\param ghost_index The index of the ghost cell
	\param time The time
	\param face_index The index of the face with the ghost point
	\param tracerstickernames The names of the tracers and stickers
	\return The gradient of the ghost cell
	*/
	virtual Slope3D GetGhostGradient(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Slope3D>& gradients, size_t ghost_index, double time, size_t face_index,
		TracerStickerNames const& tracerstickernames) const = 0;

	//! \brief Virtual destructor
	virtual ~Ghost3D(void);
	/*!
	\brief Finds the indeces of the outer faces
	\param tess The tessellation
	\return The indeces of the outer faces and whether the ghost is the first neighbor (1) or the second (2)
	*/
	vector<std::pair<size_t, size_t> > GetOuterFacesIndeces(Tessellation3D const& tess)const;
};

class RigidWallGenerator3D : public Ghost3D
{
public:
	void operator() (const Tessellation3D& tess,
		const vector<ComputationalCell3D>& cells, double time, TracerStickerNames const&
		tracerstickernames, boost::container::flat_map<size_t, ComputationalCell3D> & res) const;

	Slope3D GetGhostGradient(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Slope3D>& gradients, size_t ghost_index, double time, size_t face_index,
		TracerStickerNames const& tracerstickernames) const;
};

class FreeFlowGenerator3D : public Ghost3D
{
public:
	void operator() (const Tessellation3D& tess,
		const vector<ComputationalCell3D>& cells, double time, TracerStickerNames const&
		tracerstickernames, boost::container::flat_map<size_t, ComputationalCell3D> &res) const;

	Slope3D GetGhostGradient(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Slope3D>& gradients, size_t ghost_index, double time, size_t face_index,
		TracerStickerNames const& tracerstickernames) const;
};

class ConstantPrimitiveGenerator3D : public Ghost3D
{
private:
	const ComputationalCell3D cell_;
public:
	ConstantPrimitiveGenerator3D(ComputationalCell3D const& cell);

	void operator() (const Tessellation3D& tess,
		const vector<ComputationalCell3D>& cells, double time, TracerStickerNames const&
		tracerstickernames, boost::container::flat_map<size_t, ComputationalCell3D> &res) const;

	Slope3D GetGhostGradient(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Slope3D>& gradients, size_t ghost_index, double time, size_t face_index,
		TracerStickerNames const& tracerstickernames) const;
};

class SeveralGhostGenerator3D : public Ghost3D
{
public:
	class GhostCriteria3D
	{
	public:
		/*!
		\brief Chooses the ghost generator
		\param tess The tessellation
		\param index The index of the mesh point to calculate for
		\return The index in the ghost generator vector to choose from
		*/
		virtual size_t GhostChoose(Tessellation3D const& tess, size_t index)const = 0;

		virtual ~GhostCriteria3D(void);
	};

	/*! \brief Class constructor
	\param ghosts List of ghost generators
	\param ghostchooser Criteria for when to use each ghost generator
	*/
	SeveralGhostGenerator3D(vector<Ghost3D*> ghosts, GhostCriteria3D const& ghostchooser);

	void operator() (const Tessellation3D& tess,
		const vector<ComputationalCell3D>& cells, double time, TracerStickerNames const&
		tracerstickernames, boost::container::flat_map<size_t, ComputationalCell3D> &res) const;

	Slope3D GetGhostGradient(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Slope3D>& gradients, size_t ghost_index, double time, size_t face_index,
		TracerStickerNames const& tracerstickernames) const;
private:
	vector<Ghost3D*> ghosts_;
	GhostCriteria3D const& ghost_chooser_;
};



#endif // GHOST3D_HPP
