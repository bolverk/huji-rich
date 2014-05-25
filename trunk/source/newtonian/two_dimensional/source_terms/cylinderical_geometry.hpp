/*! \file cylinderical_geometry.hpp
  \brief Cylinderical geometry by means of source terms
  \author Almog Yalinewich
 */ 

#ifndef CYLINDERICAL_GEOMETRY_HPP
#define CYLINDERICAL_GEOMETRY_HPP 1

#include "../SourceTerm.hpp"
#include "../../../tessellation/geometry.hpp"

//! \brief Geometric source terms for cylindrical geometry
class CylindericalGeometry: public SourceTerm
{
public:

	/*! \brief Class constructor
	\param origin Axes origin
	\param direction Azimuthal direction (axis of cylinderical symmetry, aka the "z" axis)
	*/
	CylindericalGeometry(Vector2D const& origin,
		Vector2D const& direction);

	Conserved Calculate
		(Tessellation const& tess,
		vector<Primitive> const& cells,
		int point,vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const& hbc,
		vector<vector<double> > const &tracer,vector<double> &dtracer,
		vector<double> const& lengthes,double t,
		double dt);

private:
	const Vector2D origin_;
	const Vector2D direction_;
};

/*! \brief Measures distant of a point from a line
  \param point point
  \param origin Point on the line
  \param direction Slope of the line
  \return Distance
 */
double distance_from_line(Vector2D const& point,Vector2D const& origin,
			  Vector2D const& direction);

/*! \brief Cross product of a vector in the x,y plane with a unit vector in the z direction
  \param v Two dimensional vector
  \return Two dimensional vector
 */
Vector2D cross_z(Vector2D const& v);

#endif // CYLINDERICAL_GEOMETRY_HPP
