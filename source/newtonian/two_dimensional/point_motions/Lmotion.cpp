#include "Lmotion.hpp"
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

LMotion::LMotion(LinearGaussImproved const& interp, EquationOfState const& eos,EdgeVelocityCalculator const& evc,
	vector<string> skip_keys) :interp_(interp), eos_(eos),evc_(evc),skip_key_(skip_keys){}

namespace
{
	bool CheckSkip(vector<size_t> const& skip_indeces, Edge const& edge,vector<ComputationalCell> const& cells)
	{
		if (skip_indeces.empty())
			return false;
		for (size_t i = 0; i < skip_indeces.size(); ++i)
			if (cells[static_cast<size_t>(edge.neighbors.first)].stickers[skip_indeces[i]]||
				cells[static_cast<size_t>(edge.neighbors.second)].stickers[skip_indeces[i]])
				return true;
		return false;
	}
}

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
	vector<double> ws(Nedges, 0),edge_length(Nedges),density_ratios(Nedges,0);
	vector<size_t> skip_indeces;
	for (size_t i = 0; i < skip_key_.size(); ++i)
	{
		size_t loc = static_cast<size_t>(binary_find(tracerstickernames.sticker_names.begin(), tracerstickernames.sticker_names.end(), 
			skip_key_[i]) - tracerstickernames.sticker_names.begin());
		if (loc >= tracerstickernames.sticker_names.size())
			throw("Can not find skip key");
		skip_indeces.push_back(loc);
	}


#ifdef RICH_MPI
	MPI_exchange_data(tess, res, true);
#endif
	vector<Vector2D> edge_vel = evc_(tess, res);
	vector<Vector2D> normals(Nedges);
	vector<bool> to_calc(Nedges, true);
	for (size_t j = 0; j < Nedges; ++j)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(j));
		edge_length[j] = edge.GetLength();
		Primitive left = convert_to_primitive(edge_values[j].first, eos_, tracerstickernames);
		Primitive right = convert_to_primitive(edge_values[j].second, eos_, tracerstickernames);
		Vector2D p = normalize(Parallel(edge));
		normals[j] = normalize(tess.GetMeshPoint(edge.neighbors.second) -	tess.GetMeshPoint(edge.neighbors.first));
		left.Velocity = Vector2D(ScalarProd(left.Velocity, normals[j]), ScalarProd(left.Velocity, p));
		right.Velocity = Vector2D(ScalarProd(right.Velocity, normals[j]), ScalarProd(right.Velocity, p));
		ws[j] = GetWs(left,right);
		if ((tess.GetOriginalIndex(edge.neighbors.first) == tess.GetOriginalIndex(edge.neighbors.second))||
			CheckSkip(skip_indeces,edge,cells))
		{
			to_calc[j] = false;
			density_ratios[j] = 1;
		}
		else
			density_ratios[j]=cells[static_cast<size_t>(edge.neighbors.first)].density /
				cells[static_cast<size_t>(edge.neighbors.second)].density;
		density_ratios[j] = std::max(density_ratios[j], 1.0 / density_ratios[j]);
		if (edge.neighbors.first < static_cast<int>(N))
			CellLength[edge.neighbors.first] += (density_ratios[j] * edge_length[j])*Vector2D(std::abs(p.y), std::abs(p.x));
		if (edge.neighbors.second < static_cast<int>(N))
			CellLength[edge.neighbors.second] += (density_ratios[j] * edge_length[j])*Vector2D(std::abs(p.y), std::abs(p.x));
	}
	for (size_t i = 0; i < N; ++i)
		res[i] = cells[i].velocity;
	for (size_t i = 0; i < Niter; ++i)
	{
		temp.assign(N,Vector2D(0,0));
#ifdef RICH_MPI
		MPI_exchange_data(tess, res, true);
#endif
		edge_vel = evc_(tess, res);
		for (size_t j = 0; j < Nedges; ++j)
		{
			if (!to_calc[j])
				continue;
			double l = edge_length[j];
			//Vector2D p = normalize(Parallel(edge));
			double density_ratio = density_ratios[j];
			double v = ScalarProd(normals[j], edge_vel.at(j));
			double cur_ws=ws[j]-v;
			Edge const& edge = tess.GetEdge(static_cast<int>(j));
			if (edge.neighbors.first < static_cast<int>(N))
			{
				temp[edge.neighbors.first] += (density_ratio*l*cur_ws) * normals[j];
				//CellLength[edge.neighbors.first] += density_ratio*l;
				//CellLength[edge.neighbors.first] += density_ratio*l*Vector2D(std::abs(p.x),std::abs(p.y));
			}
			if (edge.neighbors.second < static_cast<int>(N))
			{
				temp[edge.neighbors.second] += (density_ratio*l*cur_ws) * normals[j];
				//CellLength[edge.neighbors.second] += density_ratio*l;
				//CellLength[edge.neighbors.second] += density_ratio*l*Vector2D(std::abs(p.x), std::abs(p.y));
			}
		}

		for (size_t j = 0; j < N; ++j)
		{
			res[j].x += (1.75 / CellLength[j].x)*temp[j].x;
			res[j].y += (1.75 / CellLength[j].y)*temp[j].y;
		}
			//res[j] += Vector2D(temp[j].x / CellLength[j].y, temp[j].y / CellLength[j].x)*abs(CellLength[j]);	
	}
	return res;
}

vector<Vector2D> LMotion::ApplyFix(Tessellation const& /*tess*/, vector<ComputationalCell> const& /*cells*/, double /*time*/,
double /*dt*/, vector<Vector2D> const& velocities, TracerStickerNames const& /*tracerstickernames*/)const
{
	return velocities;
}
