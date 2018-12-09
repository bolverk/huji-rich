#include "flux_calculator_3d.hpp"

FluxCalculator3D::~FluxCalculator3D(void) {}

namespace
{
	void AddTracers(ComputationalCell3D const& left, ComputationalCell3D const& right, Conserved3D &res)
	{
		size_t ntracers = left.tracers.size();
//		res.tracers.resize(ntracers);
		for (size_t i = 0; i < ntracers; ++i)
			res.tracers[i] = (res.mass>0 ? left.tracers[i] : right.tracers[i])*res.mass;
	}
}

void RotateSolveBack3D(Vector3D const& normal, ComputationalCell3D const& left, ComputationalCell3D const& right,
	Vector3D const& face_velocity,RiemannSolver3D const& rs, Conserved3D &res,EquationOfState const& eos,TracerStickerNames const& tsn)
{
	res = rs(left, right, ScalarProd(normal, face_velocity),eos,tsn,normal);
	// Add tracers
	AddTracers(left, right, res);
}
