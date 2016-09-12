/*! \file cylinderical_geometry.hpp
\brief Cylinderical geometry by means of source terms
\author Almog Yalinewich
*/

#ifndef CYLINDERICAL_GEOMETRY_HPP
#define CYLINDERICAL_GEOMETRY_HPP 1

#include "../SourceTerm.hpp"
#include "../../../tessellation/geometry.hpp"
#include "../../common/equation_of_state.hpp"

//! \brief Geometric source terms for cylindrical geometry
class CylindericalGeometry : public SourceTerm
{
public:

	/*! \brief Class constructor
	\param origin Axes origin
	\param direction Radial direction (a aka the "r" axis)
	\param eos THe equation of state
	*/
	CylindericalGeometry(Vector2D const& origin,Vector2D const& direction,EquationOfState const& eos);

	vector<Extensive> operator()
		(const Tessellation& tess,
			const PhysicalGeometry& pg,
			const CacheData& cd,
			const vector<ComputationalCell>& cells,
			const vector<Extensive>& fluxes,
			const vector<Vector2D>& point_velocities,
			const double t,
			TracerStickerNames const& tracerstickernames) const;


private:
	const Vector2D origin_;
	const Vector2D direction_;
	EquationOfState const& eos_;
};

/*! \brief Measures distant of a point from a line
\param point point
\param origin Point on the line
\param direction Slope of the line
\return Distance
*/
double distance_from_line(Vector2D const& point, Vector2D const& origin,
	Vector2D const& direction);

/*! \brief Cross product of a vector in the x,y plane with a unit vector in the z direction
\param v Two dimensional vector
\return Two dimensional vector
*/
Vector2D cross_z(Vector2D const& v);

#endif // CYLINDERICAL_GEOMETRY_HPP