#include <iostream>
#include <fstream>
#include "source/tessellation/calc_face_vertex_velocity.hpp"

using namespace std;

namespace {
  ostream& operator<<(ostream& os, const Vector2D& v)
  {
    os << v.x << " " << v.y;
    return os;
  }
}

int main(void)
{
  ofstream f("output.txt");

  // Translation
  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(2,3),
				 Vector2D(0,-1),Vector2D(2,3),
				 Vector2D(1,0)) << endl;

  // Stretch
  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(0,1),
				 Vector2D(0,-1),Vector2D(0,-1),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(0,-1),
				 Vector2D(0,-1),Vector2D(0,1),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,-1),Vector2D(0,-1),
				 Vector2D(0,1),Vector2D(0,1),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,-1),Vector2D(0,-1),
				 Vector2D(0,1),Vector2D(0,1),
				 Vector2D(2,0)) << endl;

  // Rotation
  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(-1,0),
				 Vector2D(0,-1),Vector2D(1,0),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(1,0),
				 Vector2D(0,-1),Vector2D(-1,0),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,-1),Vector2D(1,0),
				 Vector2D(0,1),Vector2D(-1,0),
				 Vector2D(1,0)) << endl;

  f << calc_face_vertex_velocity(Vector2D(0,1),Vector2D(-1,0),
				 Vector2D(0,-1),Vector2D(1,0),
				 Vector2D(2,0)) << endl;

  f.close();
  
  return 0;
}
