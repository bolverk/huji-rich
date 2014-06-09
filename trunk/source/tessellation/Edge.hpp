/*! \file Edge.hpp
  \brief Edge between cells
  \author Elad Steinberg
*/ 

#ifndef ___EDGE_HPP___
#define ___EDGE_HPP___

#include <vector>
#include "geometry.hpp"

using namespace std;

//! \brief Interface between two cells
class Edge
{
public:

  std::pair<Vector2D, Vector2D> vertices;
  std::pair<int, int> neighbors;

  /*! \brief Class constructor
    \param p1 Position of first vertex
    \param p2 Position of second vertex
    \param neighbor1 Index of first neighbor cell
    \param neighbor2 Index of second neighbor cell
  */
  Edge(Vector2D const& p1, Vector2D const& p2, 
       int neighbor1, int neighbor2);

  Edge(void);

  ~Edge(void);

  /*! \brief Copy constructor
    \param other Source edge
  */
  Edge(Edge const& other);

  /*! \brief Returns the length of the edge
    \return Length
  */
  double GetLength(void) const;

  //! \brief Sets the mesh generating point related to the edge. \param dim The index in the edge we want to change. \param data The index of the new mesh generating point.
  void set_friend(int dim,int data);

   /*! \brief Sets the vertex of the edge
   \param vec The new vertex
    \param index Side index (either 0 or 1)
  */
  void SetVertex(Vector2D const& vec,int index);

};

/*! \brief Exception thrown in case an invalid side index is entered
  \details The index can only be 0 or 1
*/
class InvalidSideIndexException
{
private:
  
  int _index;

public:
  /*! \brief Class constructor
    \param index The invalid input argument
  */
  InvalidSideIndexException(int index);

  /*! \brief Returns the invalid index
    \return The invalid index1
  */
  //  int GetIndex(void) const;
};

/*! \brief Calculates a unit vector parallel to an edge
  \param edge Edge
  \return Unit vector parallel to the edge
*/
Vector2D Parallel(Edge const& edge);

/*!
  \brief Calculates the distance of a point to an edge
  \param point The point
  \param edge The edge
  \returns The minimum distance between the point and the edge
*/
double DistanceToEdge(Vector2D const& point,Edge const& edge);

/*!
  \brief Calculates the intersection of two edges
  \param edge1 The first edge
  \param edge2 The second edge
  \param Intersection Gets the value of the intersection point
  \param eps The relative tolerance in units of the smaller edge
  \returns true if the intersection is on the edges and false if it is outside them
*/
bool SegmentIntersection(Edge const& edge1,Edge const& edge2,
			 Vector2D &Intersection,double eps=1e-8);

#endif	// of #ifndef ___EDGE_HPP___
