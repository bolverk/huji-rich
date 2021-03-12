/*! \file Face.hpp
  \brief Face between cells
  \author Elad Steinberg
*/

#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include <boost/container/small_vector.hpp>
#include "Vector3D.hpp"
//! \brief Container for small collection of points
typedef boost::container::small_vector<Vector3D, 10> point_vec_v;

using std::vector;

//! \brief Interface between two cells
class Face
{
public:

  //! \brief Points at the ends of the edge
  point_vec_v vertices;

  //! \brief Neighboring cells
  std::pair<std::size_t,std::size_t> neighbors;

  /*! \brief Class constructor
    \param vert Position of the vertices
    \param neighbor1 Index of first neighbor cell
    \param neighbor2 Index of second neighbor cell
  */
  Face(point_vec_v const& vert,std::size_t neighbor1,std::size_t neighbor2);

  /*! \brief Assignment operator
    \param other Source
    \return Reference to new object
   */
  Face& operator=(const Face& other);

  Face(void);

  ~Face(void);

  /*! \brief Copy constructor
    \param other Source Face
  */
  Face(Face const& other);

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
