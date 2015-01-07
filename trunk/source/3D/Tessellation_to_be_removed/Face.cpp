#include "Face.hpp"

using namespace std;
 
Face::Face(vector<Vector3D> const& vert,
	   const NeighborIndex& neighbor1,
	   const NeighborIndex& neighbor2):
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
