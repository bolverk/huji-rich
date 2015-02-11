#include "Face.hpp"
#include "../../misc/universal_error.hpp"
#include <sstream>
#include <cmath>

using namespace std;

Face::Face(vector<VectorRef> const& vert):
  vertices(vert),_neighbors() 
{
}

Face::Face(): vertices(), _neighbors() 
{
}

Face::Face(Face const& other):
  vertices(other.vertices),
  _neighbors(other._neighbors) 
{
}

double Face::GetArea(void) const
{
	double res=0;
	for(size_t i=0;i<vertices.size()-2;++i)
		res+=0.5*abs(CrossProduct(*vertices[i+1]-*vertices[0],*vertices[i+2]-*vertices[0]));
	return res;
}

bool Face::IdenticalTo(const vector<VectorRef> &otherVertices) const
{
	int size = (int)otherVertices.size();
	BOOST_ASSERT(size >= 3);

	if (vertices.size() != size)
		return false;

	// First find the index of each of our vertices in otherVertices
	vector<int> indices(size);
	for (size_t i = 0; i < size; i++)
	{
		vector<VectorRef>::const_iterator it = std::find(otherVertices.begin(), otherVertices.end(), vertices[i]);
		if (it == otherVertices.end())  // Can't find vertex in the other face
			return false;
		indices[i] = (int)std::distance(otherVertices.begin(), it);
	}

	// The indices should be either x, x+1 , x+2, ... x+size-1 (all modulu size), or
	// x, x-1, x-2, ...., x-size+1 (all modulo size) if the faces are in opposite directions.
	// So the differences between the indices should be
	// 1, 1, 1, 1, 1, -size, 1, 1, 1, 1 or -1, -1, -1, -1, size, -1, -1, -1, -1
	// with size being somewhere in there (may be in the middle, may be on one of the ends)
	bool sameDirection;
	int diff = indices[1] - indices[0];
	if (diff == 1 || diff == -size+1)
		sameDirection = true;
	else if (diff == -1 || diff == size-1)
		sameDirection = false;
	else
		return false;   // The faces are *not* identical

	bool sizeFound = false;
	for (int i = 0; i < size; i++)
	{
		int diff = indices[(i + 1) % size] - indices[i];
		if (!sameDirection)
			diff *= -1;
		if (diff != 1 && diff != -size + 1)
			return false; 
		if (diff == -size + 1)
		{
			if (sizeFound)
				return false;
			sizeFound = true;
		}
	}

	return sizeFound;
}

void Face::AddNeighbor(size_t cell, bool overlapping)
{
	// See if the neighbor already exists and overlapping hasn't changed, return gracefully
	for (auto it = _neighbors.begin(); it != _neighbors.end(); it++)
	{
		if (it->GetCell() == cell)
		{
			if (it->IsOverlapping() != overlapping)
			{
				stringstream strm;
				strm << "Can't change overlapping of neighbor " << cell << " to " << overlapping;
				throw UniversalError(strm.str());
			}
			else
				return;
		}
	}

	// If we get here, we need to add a new neighbor
	if (_neighbors.size() == 2)
	{
		stringstream strm;
		strm << "Can't add third neighrbor to cell " << cell;
		throw UniversalError(strm.str());
	}

	_neighbors.push_back(NeighborInfo(cell, overlapping));
}

Vector3D calc_centroid(const Face& face)
{
	Vector3D res;
	for(size_t i=0;i<face.vertices.size()-2;++i)
	{
		double area=0.5*abs(CrossProduct(*face.vertices[i+1]-*face.vertices[0],
			*face.vertices[i+2]-*face.vertices[0]));
		res.x+=area*(face.vertices[0]->x+face.vertices[i+2]->x+face.vertices[i+1]->x)/3.;
		res.y+=area*(face.vertices[0]->y+face.vertices[i+2]->y+face.vertices[i+1]->y)/3.;
		res.z+=area*(face.vertices[0]->z+face.vertices[i+2]->z+face.vertices[i+1]->z)/3.;
	}
	return res;
}

std::ostream& operator<<(std::ostream& output, const Face::NeighborInfo& neighbor)
{
	output << neighbor.GetCell();
	if (neighbor.IsOverlapping())
		output << "-O";
	return output;
}

static const double PI = acos(-1);  // No PI definition in the C++ standard!!!
static const double EPSILON = 1e-12;

double Face::FullAngle(const Vector3D &v1, const Vector3D &v2)
{
	double angle = CalcAngle(v1, v2);
	Vector3D cross = CrossProduct(v1, v2);

	if (cross.x < 0 ||
		cross.x == 0 && cross.y < 0 ||
		cross.x == 0 && cross.y == 0 && cross.z < 0)
	{
		// This means the angle is between Pi and 2*Pi - adjust it
		angle = 2 * PI - angle;
	}

	return angle;
}

typedef std::pair<double, VectorRef> AngledVertex;

static bool CompareAngledVertices(const AngledVertex &a1, const AngledVertex &a2)
{
	double diff = a1.first - a2.first;
	return diff < -EPSILON;
}

void Face::ReorderVertices()
{
	Vector3D center;
	// We need to sort the vectors and angles together (no standard C++ sort returns the sorting permutation)
	std::vector<std::pair<double, VectorRef>> angledVertices(vertices.size());

	for (size_t i = 0; i < vertices.size(); i++)
		center += *vertices[i];
	center = center / (int)vertices.size();

	Vector3D line0FromCenter = center - *vertices[0];
	for (size_t i = 0; i < vertices.size(); i++)
	{
		Vector3D lineFromCenter = center - *vertices[i];
		angledVertices[i].first = FullAngle(line0FromCenter, lineFromCenter);
		angledVertices[i].second = vertices[i];
	}

	std::sort(angledVertices.begin(), angledVertices.end(), CompareAngledVertices);
	
	for (size_t i = 0; i < vertices.size(); i++)  // Copy the results
		vertices[i] = angledVertices[i].second;
}
