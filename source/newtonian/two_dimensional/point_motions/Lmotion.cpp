#include "Lmotion.hpp"
#include "../simple_flux_calculator.hpp"
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

LMotion::LMotion(LinearGaussImproved const& interp, EquationOfState const& eos,EdgeVelocityCalculator const& evc) :interp_(interp), eos_(eos),
	evc_(evc){}

vector<Vector2D> LMotion::operator()(const Tessellation& tess, const vector<ComputationalCell>& cells,
double time, TracerStickerNames const& tracerstickernames) const
{
	size_t N = tess.GetPointNo();
	size_t Niter = 10;
	vector<Vector2D> res(N,Vector2D(0,0));
	vector<std::pair<ComputationalCell, ComputationalCell> > edge_values;
	edge_values.resize(static_cast<size_t>(tess.GetTotalSidesNumber()),
	pair<ComputationalCell, ComputationalCell>(cells[0], cells[0]));
	SlabSymmetry pg;
	CacheData cd(tess, pg);
	vector<Vector2D> CellLength(N);
	vector<Vector2D> temp(N);
	interp_(tess, cells, time, edge_values, tracerstickernames, cd);
	size_t Nedges=edge_values.size();
	vector<double> ws(Nedges, 0);
	//vector<double> cur_ws(ws);
#ifdef RICH_MPI
	MPI_exchange_data(tess, res, true);
#endif
	vector<Vector2D> edge_vel = evc_(tess, res);
	for (size_t j = 0; j < Nedges; ++j)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(j));
		Primitive left = convert_to_primitive(edge_values[j].first, eos_, tracerstickernames);
		Primitive right = convert_to_primitive(edge_values[j].second, eos_, tracerstickernames);
		Vector2D p = normalize(Parallel(edge));
		Vector2D n = normalize(tess.GetMeshPoint(edge.neighbors.second) -	tess.GetMeshPoint(edge.neighbors.first));
		left.Velocity = Vector2D(ScalarProd(left.Velocity, n), ScalarProd(left.Velocity, p));
		right.Velocity = Vector2D(ScalarProd(right.Velocity, n), ScalarProd(right.Velocity, p));
		ws[j] = GetWs(left,right);
	}
	for (size_t i = 0; i < N; ++i)
		res[i] = cells[i].velocity;

	for (size_t i = 0; i < Niter; ++i)
	{
		temp.assign(N,Vector2D(0,0));
		CellLength.assign(N,Vector2D(0,0));
#ifdef RICH_MPI
		MPI_exchange_data(tess, res, true);
#endif
		vector<Vector2D> edge_vel = evc_(tess, res);
		for (size_t j = 0; j < Nedges; ++j)
		{
			Edge const& edge = tess.GetEdge(static_cast<int>(j));
			double l = edge.GetLength();
			Vector2D p = normalize(Parallel(edge));
			Vector2D n = normalize(tess.GetMeshPoint(edge.neighbors.second) -	tess.GetMeshPoint(edge.neighbors.first));
			double v = ScalarProd(n, edge_vel.at(j));
			double cur_ws=ws[j]-v;
			if (edge.neighbors.first < static_cast<int>(N))
			{
				temp[edge.neighbors.first] += l*cur_ws * n;
				CellLength[edge.neighbors.first] += l*Vector2D(std::abs(p.x),std::abs(p.y));
			}
			if (edge.neighbors.second < static_cast<int>(N))
			{
				temp[edge.neighbors.second] += l*cur_ws * n;
				CellLength[edge.neighbors.second] += l*Vector2D(std::abs(p.x), std::abs(p.y));
			}
		}
 
		for (size_t j = 0; j < N; ++j)
			res[j] +=  0.4*Vector2D(temp[j].x / CellLength[j].x, temp[j].y / CellLength[j].y);
	}
	return res;
}

vector<Vector2D> LMotion::ApplyFix(Tessellation const& /*tess*/, vector<ComputationalCell> const& /*cells*/, double /*time*/,
double /*dt*/, vector<Vector2D> const& velocities, TracerStickerNames const& /*tracerstickernames*/)const
{
	return velocities;
}
