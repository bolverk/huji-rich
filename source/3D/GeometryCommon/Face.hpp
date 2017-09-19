/*! \file Face.hpp
  \brief Face between cells
  \author Elad Steinberg
*/

#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include "Vector3D.hpp"

using std::size_t;

//! \brief Interface between two cells
class Face
{
public:

  //! \brief Points at the ends of the edge
  std::vector<Vector3D> vertices;

  //! \brief Neighboring cells
  std::pair<size_t,size_t> neighbors;

  /*! \brief Class constructor
    \param vert Position of the vertices
    \param neighbor1 Index of first neighbor cell
    \param neighbor2 Index of second neighbor cell
  */
  Face(vector<Vector3D> const& vert,size_t neighbor1,size_t neighbor2);

  Face(void);

  ~Face(void);

  /*! \brief Copy constructor
    \param other Source Face
  */
  Face(Face const& other);

  Face& operator=(const Face& other);

  /*! \brief Returns the area of the face
    \return Length
  */
  double GetArea(void) const;
};

/*! \brief Calculates the centroid of aa face
  \param face The face
  \return Centroid
 */
Vector3D calc_centroid(const Face& face);

#endif	// FACE_HPP
