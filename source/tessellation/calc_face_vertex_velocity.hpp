/*! \brief Calculates the velocity of the vertices of a Voronoi face
  \author Almog Yalinewich
  \file calc_face_vertex_velocity.hpp
*/

#ifndef CALC_FACE_VERTEX_VELOCITY_HPP
#define CALC_FACE_VERTEX_VELOCITY_HPP 1

#include "geometry.hpp"

/*! \brief Calculates the velocity of the vertices of a Voronoi face
  \param p1_pos Position of the first mesh generating point
  \param p1_vel Velocity of the first mesh generating point
  \param p2_pos Position of the second mesh generating point
  \param p2_vel Velocity of the second mesh generating point
  \param vertex Vertex
  \return Velocity at vertex
 */
Vector2D calc_face_vertex_velocity(const Vector2D& p1_pos, const Vector2D& p1_vel,
				   const Vector2D& p2_pos, const Vector2D& p2_vel,
				   const Vector2D& vertex);

#endif // CALC_FACE_VERTEX_VELOCITY
