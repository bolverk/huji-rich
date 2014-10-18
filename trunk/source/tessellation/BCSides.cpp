#include "BCSides.hpp"

#include <cmath>

BCSides::BCSides(double down1, double up1,
		 double left1, double right1):
  up(up1),down(down1),left(left1),right(right1)
{
  if(down>=up||left>=right){
    WrongBCSidesOrderException eo(down, up, left, right);
    throw eo;
  }
}

BCSides::~BCSides(void)
{
}

BCSides::BCSides(void):
  up(0), down(0), left(0), right(0) {}

BCSides::BCSides(const BCSides& bc):
  up(bc.GetUp()),
  down(bc.GetDown()),
  left(bc.GetLeft()),
  right(bc.GetRight()) {}

BCSides& BCSides::operator=(const BCSides& bc)
{
	if (this == &bc)
       return *this;
	up=bc.GetUp();
	down=bc.GetDown();
	left=bc.GetLeft();
	right=bc.GetRight();
	return *this;
}

double BCSides::GetUp(void) const
{
	return up;
}

double BCSides::GetDown(void) const
{
	return down;
}

double BCSides::GetLeft(void) const
{
	return left;
}

double BCSides::GetRight(void) const
{
	return right;
}

WrongBCSidesOrderException::WrongBCSidesOrderException (double down, double up, double left, double right)
:	_up(up), _down(down), _left(left), _right(right)
{}

WrongBCSidesOrderException::~WrongBCSidesOrderException(void)
{
}

double WrongBCSidesOrderException::GetUp(void) const
{
	return _up;
}

double WrongBCSidesOrderException::GetDown(void) const
{
	return _down;
}

double WrongBCSidesOrderException::GetLeft(void) const
{
	return _left;
}

double WrongBCSidesOrderException::GetRight(void) const
{
	return _right;
}
