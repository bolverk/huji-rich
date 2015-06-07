/*! \file BCSides.hpp
  \brief Sides of the boundary of the computational domain
  \author Elad Steinberg
 */

#ifndef BCSIDES_HPP
#define BCSIDES_HPP

#include <vector>
#include "geometry.hpp"

//! \brief Positions of boundaries
class BCSides
{
private:
	 //! \brief y coordinate of upper wall
  double up;

  //! \brief y coordinate of lower wall
  double down;

  //! \brief x coordinate of left wall
  double left;

  //! \brief x coordinate of right wall
  double right;

public:
  /*! \brief Class constructor
    \param down y coordinate of the lower boundary
    \param up y coordinate of the upper boundary
    \param left x coordinate of the left boundary
    \param right x coordinate of the right boundary
   */
  BCSides(double down, double up,
	  double left, double right);

  //! \brief Default constructor, initlizes to zero.
  BCSides(void);

  //! \brief Copy constructor \param bc Other BCSides
  BCSides(const BCSides& bc);

  ~BCSides(void);

  //! \brief Assigment Operator. \param bc Other BCSides \return Copy of bc
  BCSides& operator=(const BCSides& bc);

  /*! \brief y coordinate of the upper boundary
    \return y coordinate of the upper boundary
   */
  double GetUp(void) const;

  /*! \brief y coordinate of the lower boundary
    \return y coordinate of the lower boundary
   */
  double GetDown(void) const;

  /*! \brief Returns the x coordinate of the left boundary
    \return x coordinate of the left boundary
   */
  double GetLeft(void) const;

  /*! \brief Returns the x coordinate of the right boundary
    \return x coordinate of the right boundary
   */
  double GetRight(void) const;
};

//! \brief Exception thrown if the boundaries are given in the wrong order
class WrongBCSidesOrderException
{
private:

  double _up;
  double _down;
  double _left;
  double _right;

public:

  /*! Class constructor
    \param down y coordinate of the lower boundary
    \param up y coordinate of the upper boundary
    \param left x coordinate of the left boundary
    \param right x coordinate of the right boundary
   */
  WrongBCSidesOrderException(double down, double up,
			     double left, double right);

  ~WrongBCSidesOrderException(void);

  /*! \brief Returns the y coordinate of the upper boundary
    \return y coordinate of the upper boundary
   */
  double GetUp(void) const;

  /*! \brief Returns the y coordinate of the lower boundary
    \return y coordinate of the lower boundary
   */
  double GetDown(void) const;

  /*! \brief Returns the x coordinate of the left boundary
    \return x coordinate of the left boundary
   */
  double GetLeft(void) const;

  /*! \brief Returns the x coordinate of the right boundary
    \return x coordinate of the right boundary
   */
  double GetRight(void) const;
};

#endif	// of #ifndef ___BCSIDES_HPP___
