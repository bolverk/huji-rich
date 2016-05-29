#include "LagrangianHLLC.hpp"
#include "hydrodynamic_variables.hpp"
#include "../../misc/universal_error.hpp"
#include "../../misc/utils.hpp"


Conserved LagrangianHLLC::operator()(Primitive const& left,	Primitive const& right,	double velocity) const
{
	if (is_nan(right.Velocity.x))
		throw UniversalError("Hllc::Solved entered with nan");

	const Vector2D normaldir(1, 0);
	Primitive local_left = left;
	Primitive local_right = right;

	local_left.Velocity -= velocity*normaldir;
	local_right.Velocity -= velocity*normaldir;

	double gl = local_left.Density*local_left.SoundSpeed;
	double gr = local_right.Density*local_right.SoundSpeed;
	double p = (gl*gr*(local_left.Velocity.x - local_right.Velocity.x) + gl*local_right.Pressure +
		gr*local_left.Pressure) / (gl + gr);
	double u = (gl*local_left.Velocity.x + gr*local_right.Velocity.x + (local_left.Pressure - local_right.Pressure)) / (gl + gr);
	Conserved f_gr(0,p*normaldir,p*u);

	f_gr.Energy += ScalarProd(f_gr.Momentum, velocity*normaldir);
	return f_gr;
}
