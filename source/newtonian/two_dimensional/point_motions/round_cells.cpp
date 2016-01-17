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
	Vector2D LimitShearVelocity(vector<Vector2D> &vel, Tessellation const& tess, int index)
	{
		vector<int> const& edges = tess.GetCellEdges(index);
		vector<int> neigh = tess.GetNeighbors(index);
		Vector2D const& mypoint = tess.GetMeshPoint(index);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			if (tess.GetOriginalIndex(neigh[i]) == index)
				continue;
			Edge const& edge = tess.GetEdge(edges[i]);
			const double maxdist = std::max(mypoint.distance(edge.vertices.first),
				mypoint.distance(edge.vertices.second));
			Vector2D const& otherpoint = edge.neighbors.first == index ? tess.GetMeshPoint(edge.neighbors.second) :
				tess.GetMeshPoint(edge.neighbors.first);
			Vector2D normal = otherpoint - mypoint;
			const double factor = 2 * maxdist / abs(normal);
			normal = normal/abs(normal);
			Vector2D parallel = edge.vertices.first - edge.vertices.second;
			parallel = parallel / abs(parallel);
			if (factor > 7)
			{
				return normal*ScalarProd(vel[static_cast<size_t>(index)],normal) +
				  ScalarProd((vel[static_cast<size_t>(edge.neighbors.first)]+
					      vel[static_cast<size_t>(edge.neighbors.second)])*0.5,parallel)*parallel;
			}
		}
		return vel[static_cast<size_t>(index)];
	}

	void LimitNeighborVelocity(vector<Vector2D> &vel, Tessellation const& tess,	int index, double factor)
	{
		vector<int> neigh = tess.GetNeighbors(index);
		Vector2D r = tess.GetMeshPoint(index);
		double R = tess.GetWidth(index);
		for (size_t i = 0; i<neigh.size(); ++i)
		{
			if (tess.GetOriginalIndex(neigh[i])!=index)
			{
				if (r.distance(tess.GetMeshPoint(neigh[i]))<0.1*R)
				{
					vel[static_cast<size_t>(neigh[i])] = vel[static_cast<size_t>(neigh[i])] * factor;
					return;
				}
			}
		}
	}

	void CorrectPointsOverShoot(vector<Vector2D> &v, double dt,Tessellation const& tess,OuterBoundary const&
		outer)
	{
		// check that we don't go outside grid
		size_t n = static_cast<size_t>(tess.GetPointNo());
		const double inv_dt = 1.0 / dt;
		for (size_t i = 0; i < n; ++i)
		{
			Vector2D point(tess.GetMeshPoint(static_cast<int>(i)));
			if ((v[i].x*dt * 2 + point.x)>outer.GetGridBoundary(Right))
			{
				double factor = 0.4*(outer.GetGridBoundary(Right) -	point.x)*inv_dt / abs(v[i]);
				v[i] = v[i] * factor;
				LimitNeighborVelocity(v, tess, static_cast<int>(i), factor);
			}
			if ((v[i].x*dt * 2 + point.x)<outer.GetGridBoundary(Left))
			{
				double factor = 0.4*(point.x - outer.GetGridBoundary(Left))*inv_dt / abs(v[i]);
				v[i] = v[i] * factor;
				LimitNeighborVelocity(v, tess, static_cast<int>(i), factor);
			}
			if ((v[i].y*dt * 2 + point.y)>outer.GetGridBoundary(Up))
			{
				double factor = 0.4*(outer.GetGridBoundary(Up) - point.y)*inv_dt / abs(v[i]);
				v[i] = v[i] * factor;
				LimitNeighborVelocity(v, tess, static_cast<int>(i), factor);
			}
			if ((v[i].y*dt * 2 + point.y)<outer.GetGridBoundary(Down))
			{
				double factor = 0.4*(point.y - outer.GetGridBoundary(Down))*	inv_dt / abs(v[i]);
				v[i] = v[i] * factor;
				LimitNeighborVelocity(v, tess, static_cast<int>(i), factor);
			}
		}
		return;
	}
}

Vector2D RoundCells::calc_dw(size_t i, const Tessellation& tess, const vector<ComputationalCell>& cells) const
{
	const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
	const Vector2D s = tess.GetCellCM(static_cast<int>(i));
	const double d = abs(s - r);
	const double R = tess.GetWidth(static_cast<int>(i));
	if (d < 0.9*eta_*R)
		return Vector2D(0, 0);
	const double c = std::max(eos_.dp2c(cells[i].density, cells[i].pressure,
		cells[i].tracers), abs(cells[i].velocity));
	return chi_*c*(s - r) / R;
}

Vector2D RoundCells::calc_dw(size_t i, const Tessellation& tess, double dt)const
{
	const Vector2D r = tess.GetMeshPoint(static_cast<int>(i));
	const Vector2D s = tess.GetCellCM(static_cast<int>(i));
	const double d = abs(s - r);
	const double R = tess.GetWidth(static_cast<int>(i));
	if (d < 0.9*eta_*R)
		return Vector2D(0, 0);
	const double c = d / dt;
	return chi_*c*(s - r) / R;
}

vector<Vector2D> RoundCells::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
	double time) const
{
	vector<Vector2D> res = pm_(tess, cells, time);
	if (!cold_)
	{
		for (size_t i = 0; i < res.size(); ++i)
		{
			res[i] += calc_dw(i, tess, cells);
		}
	}
	return res;
}

vector<Vector2D> RoundCells::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
	double dt, vector<Vector2D>const & velocities)const
{
	vector<Vector2D> res = pm_.ApplyFix(tess, cells, time, dt, velocities);
	res.resize(static_cast<size_t>(tess.GetPointNo()));
	if (cold_)
	{
		const size_t n = res.size();
		for (size_t i = 0; i < n; ++i)
		{
			res.at(i) += calc_dw(i, tess, dt);
		}
	}
	if (outer_.GetBoundaryType()!=Periodic)
		CorrectPointsOverShoot(res, dt, tess,outer_);
	return res;
}
