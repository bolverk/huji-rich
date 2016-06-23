#include "cylindrical_geometry.hpp"

CylindericalGeometry::CylindericalGeometry(Vector2D const& origin,Vector2D const& direction,
	EquationOfState const& eos) :
	origin_(origin), direction_(direction),eos_(eos) {}

double distance_from_line(Vector2D const& point,
	Vector2D const& origin,
	Vector2D const& direction)
{
	const double hypotenuse = abs(point - origin);
	const double side = fabs(Projection(point - origin, direction));
	return sqrt(pow(hypotenuse, 2) - pow(side, 2));
}

Vector2D cross_z(Vector2D const& v)
{
	return Vector2D(v.y, -v.x);
}

vector<Extensive> CylindericalGeometry::operator()
(const Tessellation& tess,
	const PhysicalGeometry& /*pg*/,
	const CacheData& cd,
	const vector<ComputationalCell>& cells,
	const vector<Extensive>& fluxes,
	const vector<Vector2D>& /*point_velocities*/,
	const double /*t*/,
	TracerStickerNames const& tracerstickernames) const

{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	vector<Extensive> res(N, fluxes[0]);
	for (size_t i = 0; i < N; ++i)
	{
		const double d = cells[i].density;
		const double vz = Projection(cells[i].velocity,direction_);
		const Vector2D r_hat = cross_z(direction_);
		const double vr = Projection(cells[i].velocity,r_hat);
		const double p = cells[i].pressure;
		const double e = eos_.dp2e(d, p, cells[i].tracers, tracerstickernames.tracer_names);
		const double r = distance_from_line(tess.GetCellCM(static_cast<int>(i)),origin_,direction_);
		const double r_mom = d*pow(vr, 2) / r;
		const double z_mom = d*vr*vz / r;
		const double volume = cd.volumes[i];
		res[i].mass = -volume*d*vr / r;
		res[i].momentum = -volume*Vector2D(r_mom*r_hat / abs(r_hat) + z_mom*direction_ / abs(direction_));
		res[i].energy = -volume*(vr*(p + d*(e + 0.5*pow(vr, 2) + 0.5*pow(vz, 2))) / r);
		if (!cells[i].tracers.empty())
			res[i].tracers = -d*volume*(vr / r)*cells[i].tracers;
	}
	return res;
}