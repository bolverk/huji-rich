#include "cylindrical_complementary.hpp"
#include "../../../misc/lazy_list.hpp"

namespace {
	double distance_from_axis(const Vector2D& point,
		const Axis& axis)
	{
		const double hypotenuse = abs(point - axis.origin);
		const double side = std::abs(Projection(point - axis.origin,
			axis.direction));
		return sqrt(pow(hypotenuse, 2) - pow(side, 2));
	}

	/*
	Vector2D cross_z(const Vector2D& v)
	{
	  return Vector2D(v.y,-v.x);
	}
	*/

	Vector2D remove_parallel_component(const Vector2D& v,
		const Vector2D& p)
	{
		return v - p*ScalarProd(v, p) / ScalarProd(p, p);
	}
}

CylindricalComplementary::CylindricalComplementary(const Axis& axis) :
	axis_(axis) {}

vector<Extensive> CylindricalComplementary::operator()
(const Tessellation& tess,
	const PhysicalGeometry& /*pg*/,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& /*fluxes*/,
	const vector<Vector2D>& /*point_velocities*/,
	const double /*time*/,
	TracerStickerNames const& /*tracerstickernames*/) const
{
	vector<Extensive> res(static_cast<size_t>(tess.GetPointNo()));
	const Vector2D r_hat = remove_parallel_component(tess.GetMeshPoint(0) - axis_.origin, axis_.direction);
	for (size_t i = 0; i < res.size(); ++i) 
	{
		const double p = cells[i].pressure;
		const double r = distance_from_axis
			(tess.GetCellCM(static_cast<int>(i)), axis_);
		const double volume = cd.volumes[i];
		res[i].mass = 0;
		res[i].momentum = volume*(p / r)*r_hat / abs(r_hat);
		res[i].energy = 0;
//		res[i].tracers.resize(cells[0].tracers.size(), 0);
	}
	return res;
}
