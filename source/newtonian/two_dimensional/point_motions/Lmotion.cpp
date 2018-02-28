#include "Lmotion.hpp"
#include <iostream>
#include "../simple_flux_calculator.hpp"
#include "../../../misc/utils.hpp"
#ifdef RICH_MPI
#include "../../../mpi/mpi_commands.hpp"
#endif //RICH_MPI

namespace
{
	double GetWs(Primitive const& left,Primitive const& right)
	{
		const double dl = left.Density;
		const double pl = left.Pressure;
		const double vl = left.Velocity.x;
		const double cl = left.SoundSpeed;
		const double dr = right.Density;
		const double pr = right.Pressure;
		const double vr = right.Velocity.x;
		const double cr = right.SoundSpeed;
		const double sl = std::min(vl - cl, vr - cr);
		const double sr = std::max(vl + cl, vr + cr);
		const double ss = (pr - pl + dl*vl*(sl - vl) - dr*vr*(sr - vr)) /
			(dl*(sl - vl) - dr*(sr - vr));
		return ss;
	}
}

LMotion::LMotion(SpatialReconstruction const& interp, EquationOfState const& eos,EdgeVelocityCalculator const& evc) :interp_(interp), eos_(eos),evc_(evc){}

vector<Vector2D> LMotion::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
double /*time*/, TracerStickerNames const& /*tracerstickernames*/) const
{
	size_t N = tess.GetPointNo();
	vector<Vector2D> res(N, Vector2D(0, 0));
	for (size_t i = 0; i < N; ++i)
		res[i] = cells[i].velocity;
	return res;
}

vector<Vector2D> LMotion::ApplyFix(Tessellation const& tess, vector<ComputationalCell> const& cells, double time,
double dt, vector<Vector2D> const& /*velocities*/, TracerStickerNames const& tracerstickernames)const
{
	int N = tess.GetPointNo();
	vector<Vector2D> res(N, Vector2D(0, 0));
	std::vector<double> TotalArea(N, 0);
	vector<std::pair<ComputationalCell, ComputationalCell> > edge_values;
	edge_values.resize(static_cast<size_t>(tess.GetTotalSidesNumber()),
		pair<ComputationalCell, ComputationalCell>(cells[0], cells[0]));
	SlabSymmetry pg;
	CacheData cd(tess, pg);
	interp_(tess, cells, time, edge_values, tracerstickernames, cd);
	size_t Nedges = edge_values.size();
	vector<double> ws(Nedges, 0), edge_length(Nedges);
	std::vector<Vector2D> normals(Nedges);
	for (size_t j = 0; j < Nedges; ++j)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(j));
		edge_length[j] = cd.areas[j];
		if (edge.neighbors.first < N)
			TotalArea[static_cast<size_t>(edge.neighbors.first)] += edge_length[j];
		if (edge.neighbors.second < N)
			TotalArea[static_cast<size_t>(edge.neighbors.second)] += edge_length[j];
		Primitive left = convert_to_primitive(edge_values[j].first, eos_, tracerstickernames);
		Primitive right = convert_to_primitive(edge_values[j].second, eos_, tracerstickernames);
		Vector2D p = normalize(Parallel(edge));
		normals[j] = normalize(tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first));
		left.Velocity = Vector2D(ScalarProd(left.Velocity, normals[j]), ScalarProd(left.Velocity, p));
		right.Velocity = Vector2D(ScalarProd(right.Velocity, normals[j]), ScalarProd(right.Velocity, p));
		ws[j] = GetWs(left, right);
	}
	vector<int> edges;
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	size_t Ntracers = tracerstickernames.tracer_names.size();
	for (size_t i = 0; i < static_cast<size_t>(N); ++i)
	{
		edges = tess.GetCellEdges(static_cast<int>(i));
		double A = cd.volumes[i];
		Vector2D CM = A*cd.CMs[i];
		for (size_t j = 0; j < edges.size(); ++j)
		{
			double Atemp = edge_length[edges[j]] * dt*ws[edges[j]];
			if (tess.GetEdge(edges[j]).neighbors.second == static_cast<int>(i))
				Atemp *= -1;
			Vector2D CMtemp = (tess.GetEdge(edges[j]).vertices.first + tess.GetEdge(edges[j]).vertices.second)*0.5 + ws[edges[j]] * 0.5
				*normals[edges[j]]*dt;
			CM += Atemp*CMtemp;
			A += Atemp;
		}
		res[i] = (CM / A - cd.CMs[i])/dt;
		if (indexX < Ntracers && indexY < Ntracers)
		{
			double m = 4.0*cells[i].density*cd.volumes[i]/(dt*std::sqrt(TotalArea[i]));
			Vector2D toadd(m*cells[i].tracers[indexX], m*cells[i].tracers[indexY]);
			double v = abs(res[i]);
			double t = abs(toadd);
			if (t > 0.15*v)
				 toadd = toadd*(0.15*v / t);
			res[i] += toadd;
		}
	}
	return res;
}
