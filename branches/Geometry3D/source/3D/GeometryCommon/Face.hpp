/*! \file Face.hpp
  \brief Face between cells
  \author Elad Steinberg
*/

#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include "Vector3D.hpp"
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
	std::vector<Vector3D> vertices;

	/*! \brief Class constructor
	\param vert Position of the vertices
	*/
	Face(vector<Vector3D> const& vert);
	//! \brief Default constructor
	Face();

	~Face(void);

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

#endif	// FACE_HPP
