#include "calc_face_vertex_velocity.hpp"

using namespace std;

Vector2D calc_face_vertex_velocity(const Vector2D& p1_pos, const Vector2D& p1_vel,
				   const Vector2D& p2_pos, const Vector2D& p2_vel,
				   const Vector2D& vertex)
{
  const Vector2D center = 0.5*(p2_pos + p1_pos);
  const Vector2D n = (p2_pos - p1_pos)/abs(p2_pos-p1_pos);
  const Vector2D p = zcross(p2_pos-p1_pos)/abs(zcross(p2_pos-p1_pos));

  return 0.5*(p1_vel+p2_vel)+
    p*ScalarProd(p2_vel-p1_vel,n)*ScalarProd(vertex-center,p)/abs(p2_pos-p1_pos)-
    n*ScalarProd(p2_vel-p1_vel,p)*ScalarProd(vertex-center,p)/abs(p2_pos-p1_pos);
}

