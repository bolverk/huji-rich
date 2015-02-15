/*! \file Face.hpp
  \brief Face between cells
  \author Elad Steinberg
*/

#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include "VectorRepository.hpp"
#include <boost/optional.hpp>
#include <iostream>

//! \brief Interface between two cells
class Face
{
public:
	class NeighborInfo
	{
	private:
		size_t _cell;
		bool _overlapping;

	public:
		NeighborInfo(size_t cell, bool overlapping) : _cell(cell), _overlapping(overlapping) { }

		//!\brief Returns the cell number of the neighbor
		size_t GetCell() const { return _cell;  }
		//!\brief Returns an indication whether the neighbor is 'overlapping', due to an overlapping boundary.
		bool IsOverlapping() const { return _overlapping; }
	};

private:
	std::vector<NeighborInfo> _neighbors;

public:
	//! \brief Points at the ends of the edge
	std::vector<VectorRef> vertices;

	/*! \brief Class constructor
	\param vert Position of the vertices
	*/
	Face(vector<VectorRef> const& vert);
	//Face(const vector<Vector3D> &vert);
	//! \brief Default constructor
	Face();


	/*! \brief Copy constructor
	\param other Source Face
	*/
	Face(Face const& other);

	/*! \brief Add a neighbor cell
		\param cell - the neighbor cell's index
		\param overlapping - True if the neighbor is an 'overlapping neighbor', due to overlapping boundaries.
		\remark For now this method asserts if two neighbors already exist. An exception may be added later.
		\remark Adding an already existing neighbor is legal
		\throws UniversalError { if the cell already has two neighbors }
	*/
	void AddNeighbor(size_t cell, bool overlapping=false);

	size_t NumNeighbors() const { return _neighbors.size(); }
	boost::optional<const NeighborInfo> Neighbor1() const 
	{
		if (_neighbors.size() > 0)
			return _neighbors[0];
		return boost::none;
	}
	boost::optional<const NeighborInfo> Neighbor2() const
	{
		if (_neighbors.size() > 1)
			return _neighbors[1];
		return boost::none;
	}

	boost::optional<const NeighborInfo> OtherNeighbor(int cell)
	{
		auto n1 = Neighbor1();
		if (n1.is_initialized() && n1->GetCell() == cell)
			return Neighbor2();
		return n1;
	}
	
	/*! \brief Returns the area of the face
	\return Length
	*/
	double GetArea(void) const;

	/*! \brief Reorders the vertices based on their angle from the center
	\remark This is useful when the face order wasn't constructed properly
	\remark It is assumed the vertices are all coplanar
	\remark It is assumed the vertices form a convex polygon
	*/
	void ReorderVertices();

	/*! \brief Sees if the face is identical to a list of vertices
	\param List of vertices
	\return True if this is the same face
	\remark
	    Two faces are identical if they contain the same vertices in the same order. The faces need not start
		with the same vertex, and also don't need to be the same direction, so:
		ABCD is identical to BCDA and also to CBAD but not to BACD
	*/
	bool IdenticalTo(const vector<VectorRef> &vertices) const;

private:
	/*! \brief Returns the angle between the two vectors.
	\param v1 First vector
	\param v2 Second vector
	\returns The angle between the vectors, between 0 and 2*Pi
	\remark The angle between 0 and Pi is calculated using the dot product of v1 and v2.
	  Then v1xv2 is consulted to determine if the angle is between 0 and Pi or between Pi and 2*Pi.
	  By convention, the first non-zero element of the cross product is consulted. If it's positive,
	  the angle is said to be between 0 and Pi, if it's negative it's said to be between Pi and 2*Pi
	\remark This function is used solely for ordering the faces, so this implementation is enough
	*/
	static double FullAngle(const Vector3D &v1, const Vector3D &v2);
};

/*! \brief Calculates the centroid of aa face
  \param face The face
  \return Centroid
 */
Vector3D calc_centroid(const Face& face);

//!\brief Outputs the neighbor info to the stream
std::ostream& operator<<(std::ostream& output, const Face::NeighborInfo& neighbor);

template <typename T>
std::ostream& operator<<(std::ostream& output, const boost::optional<T>& opt)
{
	if (opt)
		output << *opt;
	else
		output << "NONE";
	return output;
}

// \brief Compares to faces.
// \returns True if the faces are identical, meaning they have exactly the same vertices, regardless of ordering.
// \remarks Since faces are convex, they are the same iff they have the same vertices.
bool operator==(const Face &face1, const Face &face2);

namespace std
{
	template<>
	struct hash<Face>
	{
		typedef Face argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type &face) const
		{
			size_t value;
			std::hash<VectorRef> hasher;
			for (int i = 0; i < face.vertices.size(); i++)
				value ^= hasher(face.vertices[i]);  // Ignores the vertex ordering, as it should

			return value;
		}
	};
}

#endif	// FACE_HPP
