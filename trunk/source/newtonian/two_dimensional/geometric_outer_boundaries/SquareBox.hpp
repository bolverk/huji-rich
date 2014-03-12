/*! \file SquareBox.hpp
  \brief Square box outer boundary conditions
  \author Elad Steinberg
*/

#ifndef SQUAREBOX_HPP
#define SQUAREBOX_HPP 1

//! \brief Square box outer boundary conditions
#include "../OuterBoundary.hpp"

class SquareBox: public OuterBoundary
{
public:

  bool AreWeReflective(Edge const& edge)const;

  BoundaryType GetBoundaryType(void) const;

  double GetGridBoundary(Directions dir) const;
  /*!
    \brief Class constructor
    \param left The left coordinate
    \param right The right coordinate
    \param up The yp coordinate
    \param down The down coordinate
  */
  SquareBox(double left, double right,double up, double down);

  /*! \brief Class constructor
    \param bottom_left Coordinates of the bottom left corner
    \param top_right Coordinates of the top right corner
   */
  SquareBox(Vector2D const& bottom_left,
	    Vector2D const& top_right);

  ~SquareBox(){;}
private:
  double _left,_right,_up,_down;
};

#endif // SQUAREBOX_HPP
