/*! \file RotationForce.hpp
  \brief Class for a adding a rotation force force for cylndrical geometry assuming conservation of angular momentum
\author Elad Steinberg
*/

#ifndef ROTATIONFORCE_HPP
#define ROTATIONFORCE_HPP 1

#include "../SourceTerm.hpp"
#include "../spatial_distribution2d.hpp"
#include "../scalar_interpolation.hpp"

//! \brief Class for a adding a rotation force force for cylndrical geometry assuming conservation of angular momentum
class RotationForce : public SourceTerm
{
public:
	/*! \brief Class Constructor
	\param origin The origin of the coordinate system
	\param direction The direction of the R axis
	\param tracer_index The index of the AM tracer
	*/
	RotationForce(Vector2D const& origin,
		Vector2D const& direction,int tracer_index);

	~RotationForce(void);

	Conserved Calculate(Tessellation const* tess,
		vector<Primitive> const& cells,int point,
		vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const*hbc,
		vector<vector<double> > const &tracer,vector<double> &dtracer,
		double time,double dt);

private:
	const Vector2D origin_;
	const Vector2D direction_;
	const int tracer_index_;
};

#endif //ROTATIONFORCE_HPP
