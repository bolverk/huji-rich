#ifndef EXTERNAL_DATA_HPP
#define EXTERNAL_DATA_HPP 1

#include <string>
#include <fstream>

using std::string;
using std::ifstream;

/*! \brief Manages data from an external file in convergence runs
 */ 
class ExternalData
{
public:

  explicit ExternalData(string const& fname);

  bool getVertexMotionFlag(void) const;

  bool getInterpolationFlag(void) const;

  int getTimeIntegrationOrder(void) const;

  int getPointNumber(void) const;

private:

  bool vertex_motion_flag_;
  bool interpolation_flag_;
  int time_int_ord_;
  int points_number_;
};

#endif // EXTERNAL_DATA_HPP
