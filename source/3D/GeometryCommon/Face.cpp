#include "Face.hpp"
#include <numeric>

using std::inner_product;

using namespace std;
 
Face::Face(point_vec_v const& vert,std::size_t neighbor1,std::size_t neighbor2):
  vertices(vert),neighbors(neighbor1,neighbor2) {}

Face& Face::operator=(const Face& other)
{
  vertices = other.vertices;
  neighbors = other.neighbors;
  return *this;
}

Face::Face(void): vertices(), neighbors() {}

Face::~Face(void)
{}

Face::Face(Face const& other):
  vertices(other.vertices),
  neighbors(other.neighbors) {}

double Face::GetArea(void) const
{
  const Vector3D& ref = vertices[0];
  return inner_product(vertices.begin()+1,
		       vertices.end()-1,
		       vertices.begin()+2,
		       0,
		       [](const double x, const double y)
		       {return x+y;},
		       [&ref](const Vector3D& u, const Vector3D& v)
		       {return 0.5*fastabs(CrossProduct(u-ref,v-ref));});
}

Vector3D calc_centroid(const Face& face)
{
  const Vector3D& ref = face.vertices[0];
  return inner_product(face.vertices.begin()+1,
		       face.vertices.end()-1,
		       face.vertices.begin()+2,
		       Vector3D(0,0,0),
		       [](const Vector3D& x, const Vector3D& y)
		       {return x+y;},
		       [&ref](const Vector3D& u, const Vector3D& v)
		       {
			 const double area = 0.5*fastabs
			   (CrossProduct(u-ref, v-ref));
			 return area*(u+v+ref)/3;
		       });
}
