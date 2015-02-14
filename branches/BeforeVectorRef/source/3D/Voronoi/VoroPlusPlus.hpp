/*! \file VoroPlusPlus.hpp
\brief A Tessalation3D implementation that is a simple wrapper around the Voro++ library
\author Itay Zandbank
*/

#ifndef VOROPLUSPLUS_HPP
#define VOROPLUSPLUS_HPP

#include "../GeometryCommon/OuterBoundary3D.hpp"
#include "../GeometryCommon/Tessellation3D.hpp"
#include "TessellationBase.hpp"
#include "../GeometryCommon/Face.hpp"

#ifdef GTEST
#include "gtest/gtest.h"
#endif

class VoroPlusPlusImpl;

class VoroPlusPlus : public TessellationBase
{
public:
	VoroPlusPlus();
	virtual ~VoroPlusPlus();

	/*!
	\brief Update the tessellation
	\param points The new positions of the mesh generating points
	*/
	virtual void Update(vector<Vector3D> const& points);


	/*!
	\brief Cloning function
	*/
	virtual Tessellation3D* clone(void) const;

	/*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points (or to same cpu for single thread) ad boundary points
	\return The sent points, outer vector is the index of the outer Face and inner vector are the points sent through the face
	*/
	virtual vector<vector<size_t> >& GetDuplicatedPoints(void);
	/*!
	\brief Returns the indeces of the points that where sent to other processors as ghost points (or to same cpu for single thread) ad boundary points
	\return The sent points, outer vector is the index of the outer Face and inner vector are the points sent through the face
	*/
	virtual vector<vector<size_t> >const& GetDuplicatedPoints(void)const;
	/*!
	\brief Returns the total number of points (including ghost)
	\return The total number of points
	*/
	virtual size_t GetTotalPointNumber(void)const;

	/*!
	\brief Checks if a point is a ghost point or not
	\return True if is a ghost point, false otherwise
	*/
	virtual bool IsGhostPoint(size_t index)const;

private:
	friend class VoroPlusPlusImpl; // The PIMPL pattern: Hides the voro++ dependencies. 
								   // The class is used only in RunVoronoi, so we don't keep a pointer to it here.

	void RunVoronoi();
};

#endif // VOROPLUSPLUS_HPP