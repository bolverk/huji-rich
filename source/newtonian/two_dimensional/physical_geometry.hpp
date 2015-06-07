/*! \file physical_geometry.hpp
  \author Almog Yalinewich
  \brief Physical geometry of the grid
*/

#ifndef PHYSICAL_GEOMETRY_HPP
#define PHYSICAL_GEOMETRY_HPP 1

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <vector>
#include "../../tessellation/Edge.hpp"

using std::vector;

//! \brief Base class for physical geometry
class PhysicalGeometry
{
public:

  /*! \brief Calculates the physical area of an edge
    \param edge Cell edge
    \return Area of the edge
   */
  virtual double calcArea(const Edge& edge) const = 0;

  /*! \brief Calculates the physical volume of a cell
    \param edge_list List of edges that bound the cell
    \return volume of the cell
   */
  virtual double calcVolume(const vector<Edge>& edge_list) const = 0;

  virtual ~PhysicalGeometry(void);
};

//! \brief Slab symmetry
class SlabSymmetry: public PhysicalGeometry
{
public:

  //! \brief Class constructor
  SlabSymmetry(void);

  double calcArea(const Edge& edge) const;

  double calcVolume(const vector<Edge>& edge_list) const;
};

//! \brief Axis of revolution
class Axis
{
public:

  //! \brief Origin of the axis
  const Vector2D origin;

  //! \brief Positive direction of the axis
  const Vector2D direction;

  /*! \brief Class constructor
    \param origin_i Origin of the axis
    \param direction_i Positive direction of the axis
   */
  Axis(const Vector2D& origin_i,
       const Vector2D& direction_i);
};

//! \brief Cylindrical symmetry
class CylindricalSymmetry: public PhysicalGeometry
{
public:

  /*! \brief Class constructor
    \param origin Origin of the axis of rotation
    \param direction Positive direction of rotation axis
   */
  CylindricalSymmetry(const Vector2D& origin,
		      const Vector2D& direction);

  double calcArea(const Edge& edge) const;

  double calcVolume(const vector<Edge>& edge_list) const;

  /*! \brief Returns the axis of revolution
    \return Axis of revolution
   */
  const Axis& getAxis(void) const;

private:
  const Axis axis_;
};

#endif // PHYSICAL_GEOMETRY_HPP
