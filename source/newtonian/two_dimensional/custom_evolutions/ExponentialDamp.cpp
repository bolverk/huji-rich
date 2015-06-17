#include "ExponentialDamp.hpp"

using std::multiplies;

bool ExponentialDamp::isRelevantToInterpolation(void) const
{
	return true;
}

bool ExponentialDamp::ShouldForceTracerReset(void)const
{
	return true;
}

ExponentialDamp::ExponentialDamp(double Rmin, double Rmax, double tau, SpatialDistribution const& density,
	SpatialDistribution const& pressure, SpatialDistribution const& xvel, SpatialDistribution const& yvel,
	EquationOfState const& eos) :a_(-1.0 / (Rmax*Rmax - Rmin*Rmin)), b_(-Rmax*Rmax / (-Rmax*Rmax + Rmin*Rmin)),
	tau_(tau), density_(density), pressure_(pressure), xvel_(xvel), yvel_(yvel), eos_(eos), first_time_(true)
{}

ExponentialDamp::~ExponentialDamp(void)
{}

Conserved ExponentialDamp::CalcFlux(Tessellation const& tessellation,
	vector<Primitive> const& cells, double dt,
	SpatialReconstruction& interpolation, Edge const& edge,
	Vector2D const& facevelocity, RiemannSolver const& rs, int index,
	HydroBoundaryConditions const& bc, double time, vector<vector<double> > const& tracers)
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

Primitive ExponentialDamp::UpdatePrimitive
(vector<Conserved> const& /*conservedintensive*/,
EquationOfState const& /*eos*/,
vector<Primitive>& cells, int index, Tessellation const& tess,
double /*time*/, vector<vector<double> > const& /*tracers*/)
{
	const Vector2D point = tess.GetMeshPoint(index);
	const double r = abs(point);
	return cells[index] * (1-(a_*r*r+b_)*(1-;
}

vector<double> ExponentialDamp::UpdateTracer
(int index, vector<vector<double> >
const& /*tracers*/, vector<vector<double> > const& /*tracerchange*/, vector<Primitive> const& /*cells*/,
Tessellation const& tess, double /*time*/)
{
	return prim_.Density*tess.GetVolume(index)*tracer_;
}

vector<double> ExponentialDamp::CalcTracerFlux
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

bool ExponentialDamp::TimeStepRelevant(void)const
{
	return true;
}
