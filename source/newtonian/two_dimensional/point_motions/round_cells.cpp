
#include "round_cells.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

RoundCells::RoundCells(const PointMotion& pm, const EquationOfState& eos, OuterBoundary const& outer, double chi,
	double eta, bool cold) : pm_(pm), eos_(eos), pouter_(-1, 1, 1, -1), outer_(outer), chi_(chi), eta_(eta), cold_(cold) {}

RoundCells::RoundCells(const PointMotion& pm, const EquationOfState& eos, double chi,
	double eta, bool cold) : pm_(pm), eos_(eos), pouter_(-1, 1, 1, -1), outer_(pouter_), chi_(chi), eta_(eta),cold_(cold) {}


namespace
{
	void CorrectPointsOverShoot(vector<Vector2D> &v, double dt,Tessellation const& tess,OuterBoundary const&
		outer)
	{
		bool halfperiodic = outer.GetBoundaryType() == HalfPeriodic;
		// check that we don't go outside grid
		size_t n = static_cast<size_t>(tess.GetPointNo());
		const double inv_dt = 1.0 / dt;
		for (size_t i = 0; i < n; ++i)
		{
			Vector2D point(tess.GetMeshPoint(static_cast<int>(i)));
			double R = tess.GetWidth(static_cast<int>(i));
			if (!halfperiodic)
			{
				if ((v[i].x*dt * 2 + point.x)>outer.GetGridBoundary(Right))
				{
					double factor = 0.25*(outer.GetGridBoundary(Right) - point.x)*inv_dt / abs(v[i]);
					if (R*0.1 > (outer.GetGridBoundary(Right) - point.x))
						factor *= -0.1;
					v[i] = v[i] * factor;
				}
				if ((v[i].x*dt * 2 + point.x) < outer.GetGridBoundary(Left))
				{
					double factor = 0.25*(point.x - outer.GetGridBoundary(Left))*inv_dt / abs(v[i]);
					if (R*0.1 > (point.x - outer.GetGridBoundary(Left)))
						factor *= -0.1;
					v[i] = v[i] * factor;
				}
			}
			if ((v[i].y*dt * 2 + point.y)>outer.GetGridBoundary(Up))
			{
				double factor = 0.25*(outer.GetGridBoundary(Up) - point.y)*inv_dt / abs(v[i]);
				if (R*0.1 > (outer.GetGridBoundary(Up) - point.y))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].y*dt * 2 + point.y)<outer.GetGridBoundary(Down))
			{
				double factor = 0.25*(point.y - outer.GetGridBoundary(Down))*	inv_dt / abs(v[i]);
				if (R*0.1 > (point.y - outer.GetGridBoundary(Down)))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
		}
		return;
	}
}

Vector2D RoundCells::calc_dw(size_t i, const Tessellation& tess, const vector<ComputationalCell>& cells,
	TracerStickerNames const& tracerstickernames) const
{
	const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
	const Vector2D s = tess.GetCellCM(static_cast<int>(i));
	const double d = abs(s - r);
	const double R = tess.GetWidth(static_cast<int>(i));
	if (d < 0.9*eta_*R)
		return Vector2D(0, 0);
	const double c = std::max(eos_.dp2c(cells[i].density, cells[i].pressure,
		cells[i].tracers,tracerstickernames.tracer_names), abs(cells[i].velocity));
	return chi_*c*(s - r) / R;
}

Vector2D RoundCells::calc_dw(size_t i, const Tessellation& tess, double dt,vector<ComputationalCell> const& cells,
	TracerStickerNames const& tracerstickernames)const
{
	const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
	const Vector2D s = tess.GetCellCM(static_cast<int>(i));
	const double d = abs(s - r);
	const double R = tess.GetWidth(static_cast<int>(i));
	if (d < 0.9*eta_*R)
		return Vector2D(0, 0);
	vector<int> neigh = tess.GetNeighbors(static_cast<int>(i));
	size_t N = neigh.size();
	double cs =std::max(abs(cells[i].velocity), eos_.dp2c(cells[i].density, cells[i].pressure,
		cells[i].tracers,tracerstickernames.tracer_names));
	for (size_t j = 0; j < N; ++j)
	{
		if (tess.GetOriginalIndex(neigh[j]) == static_cast<int>(i))
			continue;
		cs = std::max(cs, eos_.dp2c(cells[static_cast<size_t>(neigh[j])].density, cells[static_cast<size_t>(neigh[j])].pressure,
			cells[static_cast<size_t>(neigh[j])].tracers));
		cs = std::max(cs, abs(cells[static_cast<size_t>(neigh[j])].velocity));
	}
	const double c_dt = d / dt;
	return chi_*std::min(c_dt,cs)*(s - r) / R;
}

vector<Vector2D> RoundCells::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
	double time, TracerStickerNames const& tracerstickernames) const
{
	vector<Vector2D> res = pm_(tess, cells, time,tracerstickernames);
	if (!cold_)
	{
		for (size_t i = 0; i < res.size(); ++i)
		{
			res[i] += calc_dw(i, tess, cells,tracerstickernames);
		}
	}
	return res;
}

vector<Vector2D> RoundCells::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	double dt, vector<Vector2D>const & velocities, TracerStickerNames const& tracerstickernames)const
{
	vector<Vector2D> res = pm_.ApplyFix(tess, cells, time, dt, velocities,tracerstickernames);
	res.resize(static_cast<size_t>(tess.GetPointNo()));
	if (cold_)
	{
		const size_t n = res.size();
		for (size_t i = 0; i < n; ++i)
		{
			res.at(i) += calc_dw(i, tess, dt,cells,tracerstickernames);
		}
	}
	if (outer_.GetBoundaryType()!=Periodic)
		CorrectPointsOverShoot(res, dt, tess,outer_);
	return res;
}
