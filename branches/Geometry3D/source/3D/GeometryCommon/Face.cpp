#include "Face.hpp"
#include <cassert>

using namespace std;

Face::Face(vector<Vector3D> const& vert,size_t neighbor1,size_t neighbor2):
  vertices(vert),neighbors(neighbor1,neighbor2) {}

Face::Face(void): vertices(), neighbors() {}

  Face::~Face(void)
  {}

Face::Face(Face const& other):
  vertices(other.vertices),
  neighbors(other.neighbors) {}

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

bool Face::AddNeighbor(size_t cell)
{
	if (neighbors.first == cell || neighbors.second == cell)
		return true;
	assert(neighbors.first == NO_NEIGHBOR || neighbors.second == NO_NEIGHBOR);
	if (neighbors.first == NO_NEIGHBOR)
		neighbors.first = cell;
	else if (neighbors.second == NO_NEIGHBOR)
		neighbors.second = cell;
	else
		return false;

	return true;
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
