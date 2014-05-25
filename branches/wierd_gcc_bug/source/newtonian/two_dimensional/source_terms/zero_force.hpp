/*! \brief Zero external force module
\author Elad Steinberg
*/
#ifndef ZEROFORCE_HPP
#define ZEROFORCE_HPP 1

#include "../SourceTerm.hpp"

class ZeroForce: public SourceTerm
{
public:
	Conserved Calculate(Tessellation const* tess,
		vector<Primitive> const& cells,int point,
		vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const *hbc,
		vector<vector<double> > const &tracer_extensive,vector<double> &dtracer,
		double t,
		double dt);
};
#endif // ZEROFORCE_HPP
