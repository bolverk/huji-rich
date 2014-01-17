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
		(Tessellation const* tess,
		vector<Primitive> const& cells,
		int point,vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const*hbc,
		vector<vector<double> > const &tracer,vector<double> &dtracer,
		double t,
		double dt);

private:
	const Vector2D origin_;
	const Vector2D direction_;
};

double distance_from_line(Vector2D const& point,Vector2D const& origin,
	Vector2D const& direction);

Vector2D cross_z(Vector2D const& v);

#endif // CYLINDERICAL_GEOMETRY_HPP
