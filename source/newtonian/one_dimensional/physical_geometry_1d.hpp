/*! \file physical_geometry_1d.hpp
  \author Almog Yalinewich
  \brief Physical geometry of the grid
*/

#ifndef PHYSICAL_GEOMETRY_1D
#define PHYSICAL_GEOMETRY_1D

//! \brief Base class for physical geometry
class PhysicalGeometry1D
{
public:

  /*! \brief Calculates the area
    \param radius Radius
    \return Area
   */
  virtual double calcArea(double radius) const = 0;

  /*! \brief Calculates the volume
    \param radius Radius
    \return Volume
   */
  virtual double calcVolume(double radius) const = 0;

  virtual ~PhysicalGeometry1D(void);
};

//! \brief Planar geometry
class SlabSymmetry1D: public PhysicalGeometry1D
{
public:

  //! \brief Class constructor
  SlabSymmetry1D(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

//! \brief Cylindrical geometry
class CylindricalSymmetry1D: public PhysicalGeometry1D
{
public:

  //! \brief Class constructor
  CylindricalSymmetry1D(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

//! \brief Spherical geometry
class SphericalSymmetry1D: public PhysicalGeometry1D
{
public:
  //! \brief Class constructor
  SphericalSymmetry1D(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

#endif // PHYSICAL_GEOMETRY_1D
