//\file CellCalculations.cpp
//\brief Performs various calculations on cells

#include "CellCalculations.hpp"
#include <vector>
#include <unordered_set>
#include "VectorRepository.hpp"

using namespace std;

//\brief Splits a cell into tetrahedra, all touching the center of the cell
std::vector<Tetrahedron> SplitCell(const std::vector<const Face *> &cell)
{
	std::vector<Tetrahedron> tetrahedra;
	Vector3D center;
	std::unordered_set<VectorRef> considered;

	// Find the center of the cell (an average of all the vertices)
	for (std::vector<const Face *>::const_iterator itFace = cell.begin(); itFace != cell.end(); itFace++)
	{
		const Face *face = *itFace;
		for (std::vector<VectorRef>::const_iterator itVertex = face->vertices.cbegin(); itVertex != face->vertices.cend(); itVertex++)
		{
			if (considered.find(*itVertex) != considered.end())  // See if we've used this vertex before
				continue;
			considered.insert(*itVertex);
			center += **itVertex;
		}
	}
	center = center / (double)considered.size();   // Average
	VectorRef centerRef(center);

	// Now create the tetrahedra, from the center to each of the faces
	for (std::vector<const Face *>::const_iterator itFace = cell.begin(); itFace != cell.end(); itFace++)
	{
		const Face *face = *itFace;
		// Split the face into trianges (face[0], face[1], face[2]), (face[0], face[2], face[3]) and so on until (face[0], face[n-2], face[n-1])
		// add center to each triangle, providing the tetrahedron
		for (size_t i = 1; i < face->vertices.size() - 1; i++)
			tetrahedra.push_back(Tetrahedron(centerRef, face->vertices[0], face->vertices[i], face->vertices[i + 1]));
	}

	return tetrahedra;
}

//\brief Calculates the volume and center-of-mass of a cell
void CalculateCellDimensions(const std::vector<const Face *> &cell, double &volume, Vector3D &centerOfMass)
{
	std::vector<Tetrahedron> tetrahedra = SplitCell(cell);
	volume = 0;
	centerOfMass = Vector3D();

	for (size_t j = 0; j < tetrahedra.size(); j++)
		volume += tetrahedra[j].volume();

	for (size_t j = 0; j < tetrahedra.size(); j++)
	{
		Vector3D weightedCenter = *tetrahedra[j].center() * tetrahedra[j].volume() / volume;
		centerOfMass += weightedCenter;
	}
}
