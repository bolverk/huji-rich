#include "Face.hpp"
#include "../../misc/universal_error.hpp"
#include <sstream>
#include <cmath>

using namespace std;

Face::Face(vector<Vector3D> const& vert):
  vertices(vert),_neighbors() 
{
}

Face::Face(void): vertices(), _neighbors() {}

  Face::~Face(void)
  {}

Face::Face(Face const& other):
  vertices(other.vertices),
  _neighbors(other._neighbors) {}

double Face::GetArea(void) const
{
	double res=0;
	for(size_t i=0;i<vertices.size()-2;++i)
		res+=0.5*abs(CrossProduct(vertices[i+1]-vertices[0],vertices[i+2]-vertices[0]));
	return res;
}

bool Face::IdenticalTo(const vector<Vector3D> otherVertices) const
{
	if (vertices.size() != otherVertices.size())
		return false;

	// This implemenation is O(N^2), but since faces contain just a few vertices, it should be fine.
	// To be more efficient (O(NlogN)) we need to be able to compare vectors.
	// Using a hash is problematic, as doubles make lousy hash keys, and the performance gain to O(N) is most
	// certainly not worh it.

	for (size_t i = 0; i < vertices.size(); i++)
	{
		bool found = false;
		for (size_t j = 0; j < otherVertices.size(); j++)
			if (vertices[i] == otherVertices[j])
			{
				found = true;
				break;
			}
		if (!found)
			return false;
	}

	return true;
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
		double area=0.5*abs(CrossProduct(face.vertices[i+1]-face.vertices[0],
			face.vertices[i+2]-face.vertices[0]));
		res.x+=area*(face.vertices[0].x+face.vertices[i+2].x+face.vertices[i+1].x)/3.;
		res.y+=area*(face.vertices[0].y+face.vertices[i+2].y+face.vertices[i+1].y)/3.;
		res.z+=area*(face.vertices[0].z+face.vertices[i+2].z+face.vertices[i+1].z)/3.;
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

typedef std::pair<double, Vector3D> AngledVertex;

static int CompareAngledVertices(const AngledVertex &a1, const AngledVertex &a2)
{
	double diff = a1.first - a2.first;
	if (diff < -EPSILON)
		return -1;
	else if (diff > EPSILON)
		return 1;
	return 0;
}

void Face::ReorderVertices()
{
	Vector3D center;
	// We need to sort the vectors and angles together (no standard C++ sort returns the sorting permutation)
	std::vector<std::pair<double, Vector3D>> angledVertices(vertices.size());

	for (size_t i = 0; i < vertices.size(); i++)
		center += vertices[i];
	center = center / vertices.size();

	Vector3D line0FromCenter = center - vertices[0];
	for (size_t i = 0; i < vertices.size(); i++)
	{
		Vector3D lineFromCenter = center - vertices[i];
		angledVertices[i].first = FullAngle(line0FromCenter, lineFromCenter);
		angledVertices[i].second = vertices[i];
	}

	std::sort(angledVertices.begin(), angledVertices.end(), CompareAngledVertices);
	
	for (size_t i = 0; i < vertices.size(); i++)  // Copy the results
		vertices[i] = angledVertices[i].second;
}
