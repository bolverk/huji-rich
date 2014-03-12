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
private:
  Vector2D _p1;
  Vector2D _p2;
  int _neighbor1;
  int _neighbor2;
public:

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

  /*! \brief Returns the vertex of the edge
    \param index Side index (either 0 or 1)
    \return Position of vertex
  */
  Vector2D GetVertex(int index) const;

  /*! \brief Returns the length of the edge
    \return Length
  */
  double GetLength(void) const;

  /*! \brief Returns the index of neighbor cells
    \param index side index (either 0 or 1)
    \return Cell index
  */
  int GetNeighbor(int index) const;

  /*! \brief Returns the x-coordinate from one of the ends of the edge. 
    \param index Which end of the edge to return. 
    \returns The chosen x-coordinate.
  */
  double get_x(int index) const;

  /*! \brief Returns the y-coordinate from one of the ends of the edge. 
    \param index Which end of the edge to return. 
    \returns The chosen y-coordinate.
  */
  double get_y(int index)const;

  //! \brief Sets the x-coordinate of a chosen end of the edge. \param point The end to change. \param data The new x-coordinate.
  void set_x(int point,double data);

  //! \brief Sets the y-coordinate of a chosen end of the edge. \param point The end to change. \param data The new y-coordinate.
  void set_y(int point,double data);

  //! \brief Sets the mesh generating point related to the edge. \param dim The index in the edge we want to change. \param data The index of the new mesh generating point.
  void set_friend(int dim,int data);

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
  \returns true if the intersection is on the edges and false if it is outside them
*/
bool SegmentIntersection(Edge const& edge1,Edge const& edge2,
			 Vector2D &Intersection);

#endif	// of #ifndef ___EDGE_HPP___
