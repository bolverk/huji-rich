#include "ConstantPrimitiveKepler.hpp"

using std::multiplies;

bool ConstantPrimitiveKepler::isRelevantToInterpolation(void) const
{
	return true;
}

bool ConstantPrimitiveKepler::ShouldForceTracerReset(void)const
{
	return true;
}

ConstantPrimitiveKepler::ConstantPrimitiveKepler(double Mass, double density, double pressure) :
M_(Mass), rho_(density), p_(pressure)
{}

ConstantPrimitiveKepler::~ConstantPrimitiveKepler(void)
{}

Conserved ConstantPrimitiveKepler::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells, double dt,
	SpatialReconstruction& interpolation, Edge const& edge,
	Vector2D const& facevelocity, RiemannSolver const& rs, int /*index*/,
	HydroBoundaryConditions const& bc, double time, vector<vector<double> > const& /*tracers*/)
{
	if (bc.IsBoundary(edge, tessellation))
		return bc.CalcFlux(tessellation, cells, facevelocity, edge, interpolation, dt, time);
	else
	{
		Vector2D normaldir = tessellation.GetMeshPoint(edge.neighbors.second) -
			tessellation.GetMeshPoint(edge.neighbors.first);

		Vector2D paraldir = edge.vertices.second - edge.vertices.first;

		Primitive left = interpolation.Interpolate
			(tessellation, cells, dt, edge, 0, InBulk, facevelocity);
		Primitive right = interpolation.Interpolate
			(tessellation, cells, dt, edge, 1, InBulk, facevelocity);
		Conserved res(FluxInBulk(normaldir, paraldir, left, right, facevelocity, rs));

		return res;
	}
}

Primitive ConstantPrimitiveKepler::UpdatePrimitive
(vector<Conserved> const& /*conservedintensive*/,
EquationOfState const& eos,
vector<Primitive>& /*cells*/, int index, Tessellation const& tess,
double /*time*/, vector<vector<double> > const& /*tracers*/)
{
	Primitive res;
	Vector2D const& point = tess.GetCellCM(index);
	double r = abs(point);
	res.Density = rho_;
	res.Pressure = p_;
	res.Velocity = sqrt(M_*pow(r, -3.0))*Vector2D(-point.y, point.x);
	res.Energy = eos.dp2e(rho_, p_);
	res.SoundSpeed = eos.dp2c(rho_, p_);
	return res;
}

vector<double> ConstantPrimitiveKepler::UpdateTracer
(int index, vector<vector<double> >
const& tracers, vector<vector<double> > const& /*tracerchange*/, vector<Primitive> const& cells,
Tessellation const& tess, double /*time*/)
{
	return tess.GetVolume(index)*cells[index].Density*tracers[index];
}

vector<double> ConstantPrimitiveKepler::CalcTracerFlux
(Tessellation const& tess,
vector<Primitive> const& cells, vector<vector<double> > const& tracers,
double dm, Edge const& edge, int /*index*/, double dt, double /*time*/,
SpatialReconstruction const& interp, Vector2D const& vface)
{
	vector<double> res(tracers[0].size());
	if (dm>0)
	{
		res = interp.interpolateTracers(tess, cells, tracers, dt, edge, 0,
			InBulk, vface);
		transform(res.begin(), res.end(), res.begin(),
			bind1st(multiplies<double>(), dm*dt*edge.GetLength()));
	}
	else
	{
		res = interp.interpolateTracers(tess, cells, tracers, dt, edge, 1,
			Boundary, vface);
		transform(res.begin(), res.end(), res.begin(),
			bind1st(multiplies<double>(), dm*dt*edge.GetLength()));
	}
	return res;
}

bool ConstantPrimitiveKepler::TimeStepRelevant(void)const
{
	return false;
}
