/*! \file SquareBox.hpp
  \brief Square box outer boundary conditions
  \author Elad Steinberg
*/

#ifndef SQUAREBOX_HPP
#define SQUAREBOX_HPP 1

//! \brief Square box outer boundary conditions
#include "../OuterBoundary.hpp"

using std::pair;

//! \brief Square frame for the tessellation
class SquareBox: public OuterBoundary
{
public:

  bool AreWeReflective(Edge const& edge)const override;

  BoundaryType GetBoundaryType(void) const override;

  bool PointIsReflective(Vector2D const& point)const override;

  double GetGridBoundary(Directions dir) const override;
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

  /*! \brief Returns the coordinates of the lower left and top right of the square frame
    \return pair of vectors: first is lower left and second is top right
   */
  pair<Vector2D,Vector2D> getBoundary(void) const;

  ~SquareBox(void) override;
private:
  double _left,_right,_up,_down;
};

#endif // SQUAREBOX_HPP
