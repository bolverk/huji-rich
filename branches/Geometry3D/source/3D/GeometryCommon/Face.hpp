/*! \file Face.hpp
  \brief Face between cells
  \author Elad Steinberg
*/

#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include "Vector3D.hpp"

//! \brief Interface between two cells
class Face
{
public:

	//! \brief Points at the ends of the edge
	std::vector<Vector3D> vertices;

	//! \brief Neighboring cells
	std::pair<size_t,size_t> neighbors;

	const static size_t NO_NEIGHBOR = (size_t)-1;

	/*! \brief Class constructor
	\param vert Position of the vertices
	\param neighbor1 Index of first neighbor cell
	\param neighbor2 Index of second neighbor cell
	*/
	Face(vector<Vector3D> const& vert,size_t neighbor1 = NO_NEIGHBOR,size_t neighbor2 = NO_NEIGHBOR);
	Face();

	~Face(void);

	/*! \brief Copy constructor
	\param other Source Face
	*/
	Face(Face const& other);

	/*! \brief Add a neighbor cell
		\param cell - the neighbor cell's index
		\returns true if the neighbor was set, false it it couldn't be set (this will need to change into an exception)
		\remark For now this method asserts if two neighbors already exist. An exception may be added later.
		\remark Adding an already existing neighbor is legal
	*/
	bool AddNeighbor(size_t cell);

	/*! \brief Returns the area of the face
	\return Length
	*/
	double GetArea(void) const;


	/*! \brief Sees if the face is identical to a list of vertices
	\param List of vertices
	\return True if this is the same face
	\remark
		We just make sure the list of vertices is identical to our list. Since we assume faces are
		convex, as long as all the vertices appear in both, it's the same face.
	*/
	bool IdenticalTo(const vector<Vector3D> vertices) const;
};

/*! \brief Calculates the centroid of aa face
  \param face The face
  \return Centroid
 */
Vector3D calc_centroid(const Face& face);

#endif	// FACE_HPP
