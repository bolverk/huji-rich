#include "flux_calculator_3d.hpp"

FluxCalculator3D::~FluxCalculator3D(void) {}

namespace
{
	Primitive convert_to_primitive(ComputationalCell3D const& cell,	const EquationOfState& eos,TracerStickerNames const& tsn,
		Vector3D const& normal)
	{
		const double energy = eos.dp2e(cell.density, cell.pressure, cell.tracers, tsn.tracer_names);
		const double sound_speed = eos.dp2c(cell.density, cell.pressure, cell.tracers, tsn.tracer_names);
		return Primitive(cell.density, cell.pressure, Vector2D(ScalarProd(cell.velocity, normal), 0), energy,
			sound_speed);
	}

	void AddParallelComponents(Conserved3D &res, ComputationalCell3D const&left, ComputationalCell3D const& right,
		Vector3D const& normal)
	{
		Vector3D par_velocity;
		if (res.mass > 0)
			par_velocity = left.velocity - normal*ScalarProd(left.velocity, normal);
		else
			par_velocity = right.velocity - normal*ScalarProd(right.velocity, normal);
		res.momentum += res.mass*par_velocity;
		res.energy += 0.5*res.mass*ScalarProd(par_velocity, par_velocity);
	}

	void AddTracers(ComputationalCell3D const& left, ComputationalCell3D const& right, Conserved3D &res)
	{
		size_t ntracers = left.tracers.size();
		res.tracers.resize(ntracers);
		for (size_t i = 0; i < ntracers; ++i)
			res.tracers[i] = (res.mass>0 ? left.tracers[i] : right.tracers[i])*res.mass;
	}
}

void RotateSolveBack3D(Vector3D const& normal, ComputationalCell3D const& left, ComputationalCell3D const& right,
	Vector3D const& face_velocity,RiemannSolver const& rs, Conserved3D &res,EquationOfState const& eos,TracerStickerNames const& tsn)
{
	Primitive primitive_left = convert_to_primitive(left, eos, tsn, normal);
	Primitive primitive_right = convert_to_primitive(right, eos, tsn, normal);
	Conserved res_temp = rs(primitive_left, primitive_right,ScalarProd(normal,face_velocity));
	res.mass = res_temp.Mass;
	res.energy = res_temp.Energy;
	res.momentum = res_temp.Momentum.x*normal;
	// Add parallel momentum and energy
	AddParallelComponents(res, left, right, normal);
	// Add tracers
	AddTracers(left, right, res);
}
