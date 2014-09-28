/*! \file physical_geometry.hpp
  \author Almog Yalinewich
  \brief Physical geometry of the grid
*/

#ifndef PHYSICAL_GEOMETRY_HPP
#define PHYSICAL_GEOMETRY_HPP 1

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "../../tessellation/Edge.hpp"

using std::vector;

class PhysicalGeometry
{
public:

  virtual double calcArea(const Edge& edge) const = 0;

  virtual double calcVolume(const vector<Edge>& edge_list) const = 0;

  virtual ~PhysicalGeometry(void);
};

class SlabSymmetry: public PhysicalGeometry
{
public:

  SlabSymmetry(void);

  double calcArea(const Edge& edge) const;

  double calcVolume(const vector<Edge>& edge_list) const;
};

class Axis
{
public:

  const Vector2D origin;
  const Vector2D direction;

  Axis(const Vector2D& origin_i,
       const Vector2D& direction_i);
};

class CylindricalSymmetry: public PhysicalGeometry
{
public:

  CylindricalSymmetry(const Vector2D& origin,
		      const Vector2D& direction);

  double calcArea(const Edge& edge) const;

  double calcVolume(const vector<Edge>& edge_list) const;

  const Axis& getAxis(void) const;

private:
  const Axis axis_;
};

#endif // PHYSICAL_GEOMETRY_HPP
