#include <iostream>
#include "external_data.hpp"
#include "../misc/universal_error.hpp"

ExternalData::ExternalData(string const& fname):
  vertex_motion_flag_(false),
  interpolation_flag_(false),
  time_int_ord_(0),
  points_number_(0)
{
  char buf_c;
  int buf_i;

  ifstream f(fname.c_str());
  f >> buf_c;
  if('e'==buf_c)
    vertex_motion_flag_ = false;
  else if('l'==buf_c)
    vertex_motion_flag_ = true;
  else
    throw UniversalError(string("First argument in external file should either be e or l (short for eulerian or lagrangian) instead, it is ")+buf_c);

  f >> buf_c;
  if('f'==buf_c)
    interpolation_flag_ = false;
  else if('s'==buf_c)
    interpolation_flag_ = true;
  else
    throw UniversalError("Second argument in external file should either be f or s (short for first or second order");

  f >> buf_i;
  time_int_ord_ = buf_i;

  f >> buf_i;
  if(buf_i>=10)
    points_number_ = buf_i;
  else
    throw UniversalError("Number of point must be larger than 9");

  f.close();
}

bool ExternalData::getVertexMotionFlag
(void) const
{
  return vertex_motion_flag_;
}

bool ExternalData::getInterpolationFlag
(void) const
{
  return interpolation_flag_;
}

int ExternalData::getTimeIntegrationOrder
(void) const
{
  return time_int_ord_;
}

int ExternalData::getPointNumber(void) const
{
  return points_number_;
}
