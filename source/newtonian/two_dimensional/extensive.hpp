#ifndef EXTENSIVE_HPP
#define EXTENSIVE_HPP 1

#include <map>
#include "../../tessellation/geometry.hpp"

class Extensive
{
public:
  
  double mass;
  double energy;
  Vector2D momentum;

  std::map<std::string,double> tracers;

  Extensive& operator-=(const Extensive& diff);

  Extensive& operator+=(const Extensive& diff);

  Extensive(void);
};

Extensive operator*(const double s,
		    const Extensive& e);

#endif // EXTENSIVE_HPP
