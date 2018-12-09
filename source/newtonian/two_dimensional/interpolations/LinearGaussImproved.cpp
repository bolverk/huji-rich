#include "LinearGaussImproved.hpp"
#include "../../../misc/utils.hpp"
#ifdef RICH_MPI
#include "../../../mpi/mpi_commands.hpp"
#endif //RICH_MPI



namespace
{
	void GetNeighborMesh(Tessellation const& tess, vector<const Edge *> const& edges, size_t cell_index,
		vector<Vector2D> &res)
	{
		res.resize(edges.size());
		const size_t nloop = edges.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			const int neigh0 = edges[i]->neighbors.first;
			const int neigh1 = edges[i]->neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
				res[i] = tess.GetMeshPoint(neigh1);
			else
				res[i] = tess.GetMeshPoint(neigh0);
		}
	}

	void GetNeighborCM(Tessellation const& tess, vector<const Edge*> const& edges, size_t cell_index,
		vector<Vector2D> &res)
	{
		res.resize(edges.size());
		const size_t nloop = edges.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			const int neigh0 = edges[i]->neighbors.first;
			const int neigh1 = edges[i]->neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
				res[i] = tess.GetCellCM(neigh1);
			else
				res[i] = tess.GetCellCM(neigh0);
		}
	}

	vector<ComputationalCell const*> GetNeighborCells(vector<const Edge *> const& edges, size_t cell_index,
		vector<ComputationalCell> const& cells)
	{
		vector<ComputationalCell const*> res(edges.size());
		const size_t nloop = edges.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			size_t other_cell = (edges[i]->neighbors.first == static_cast<int>(cell_index)) ? static_cast<size_t>
				(edges[i]->neighbors.second) : static_cast<size_t> (edges[i]->neighbors.first);
			res[i] = &cells.at(other_cell);
		}
		return res;
	}

	void GetEdgeList(Tessellation const& tess,
		vector<int> const& edge_indices, vector<const Edge*> &res)
	{
		res.clear();
		res.reserve(edge_indices.size());
		const size_t nloop = edge_indices.size();
		for (size_t i = 0; i < nloop; ++i)
			res.push_back(&tess.GetEdge(edge_indices[i]));
	}

	void calc_naive_slope(ComputationalCell const& cell,
		Vector2D const& center, Vector2D const& cell_cm, double cell_volume, vector<ComputationalCell const*> const& neighbors,
		vector<Vector2D> const& neighbor_centers, vector<Vector2D> const& neigh_cm, vector<const Edge*> const& edge_list,
		Slope &res,Slope &vec_compare)
	{
		size_t n = edge_list.size();
		if (n > 20)
		{
			UniversalError eo("Cell has too many neighbors");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			throw eo;
		}
		// Create the matrix to invert and the vector to compare
		vector<double> m(4, 0);
		for (size_t i = 0; i < edge_list.size(); ++i)
		{
			const Vector2D c_ij = CalcCentroid(*edge_list[i]) - 0.5*(neigh_cm[i] + cell_cm);
			const double e_length = edge_list[i]->GetLength();
			const Vector2D r_ij = normalize(neighbor_centers[i] - center)*e_length;
			m[0] -= c_ij.x*r_ij.x;
			m[1] -= c_ij.y*r_ij.x;
			m[2] -= c_ij.x*r_ij.y;
			m[3] -= c_ij.y*r_ij.y;
			if (i == 0)
			{
				ReplaceComputationalCell(vec_compare.xderivative, cell);
				ReplaceComputationalCell(vec_compare.yderivative, cell);
				vec_compare.xderivative *= r_ij.x*0.5;
				vec_compare.yderivative *= r_ij.y*0.5;
			}
			else
			{
				ComputationalCellAddMult(vec_compare.xderivative, cell, r_ij.x*0.5);
				ComputationalCellAddMult(vec_compare.yderivative, cell, r_ij.y*0.5);
			}
			ComputationalCellAddMult(vec_compare.yderivative, *neighbors[i], r_ij.y*0.5);
			ComputationalCellAddMult(vec_compare.xderivative, *neighbors[i], r_ij.x*0.5);
		}
		m[0] += cell_volume;
		m[3] += cell_volume;
		// Find the det
		const double det = m[0] * m[3] - m[1] * m[2];
		// Check none singular
		if (std::abs(det) < 1e-10*cell_volume*cell_volume)
		{
			UniversalError eo("Singular matrix");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			eo.AddEntry("Cell volume", cell_volume);
			eo.AddEntry("Det was", det);
			throw eo;
		}
		// Invert the matrix
		vector<double> m_inv(4);
		const double det_inv = 1.0 / det;
		m_inv[0] = m[3] * det_inv;
		m_inv[1] = -m[1] * det_inv;
		m_inv[2] = -m[2] * det_inv;
		m_inv[3] = m[0] * det_inv;
		// Calculate the gradient
		ReplaceComputationalCell(res.xderivative, vec_compare.xderivative);
		ReplaceComputationalCell(res.yderivative, vec_compare.yderivative);
		res.xderivative *= m_inv[0];
		res.yderivative *= m_inv[3];
		ComputationalCellAddMult(res.xderivative, vec_compare.yderivative, m_inv[1]);
		ComputationalCellAddMult(res.yderivative, vec_compare.xderivative, m_inv[2]);
	}

	double PressureRatio(ComputationalCell cell, vector<ComputationalCell const*> const& neigh)
	{
		double res = 1;
		double p = cell.pressure;
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (p > neigh[i]->pressure)
				res = std::min(res, neigh[i]->pressure / p);
			else
				res = std::min(res, p / neigh[i]->pressure);
		}
		return res;
	}

	bool is_shock(Slope const& naive_slope, double cell_width, double shock_ratio,
		ComputationalCell const& cell, vector<ComputationalCell const*> const& neighbor_list, double pressure_ratio, double cs)
	{
		const bool cond1 = (naive_slope.xderivative.velocity.x + naive_slope.yderivative.velocity.y)*
			cell_width < (-shock_ratio)*cs;
		const bool cond2 = PressureRatio(cell, neighbor_list) < pressure_ratio;
		return cond1 || cond2;
	}

	ComputationalCell interp(ComputationalCell const& cell, Slope const& slope,
		Vector2D const& target, Vector2D const& cm)
	{
		ComputationalCell res(cell);
		ComputationalCellAddMult(res, slope.xderivative, target.x - cm.x);
		ComputationalCellAddMult(res, slope.yderivative, target.y - cm.y);
		return res;
	}

	void interp2(ComputationalCell &res, Slope const& slope,
		Vector2D const& target, Vector2D const& cm)
	{
		ComputationalCellAddMult(res, slope.xderivative, target.x - cm.x);
		ComputationalCellAddMult(res, slope.yderivative, target.y - cm.y);
	}

	void slope_limit(ComputationalCell const& cell, Vector2D const& cm,
		vector<ComputationalCell const*> const& neighbors, vector<const Edge*> const& edge_list,
		Slope &slope,
		ComputationalCell &cmax,
		ComputationalCell &cmin,
		ComputationalCell &maxdiff,
		ComputationalCell &mindiff,
		TracerStickerNames const& tracerstickernames,
		string const& skip_key,
		Tessellation const& tess)
	{
		ReplaceComputationalCell(cmax, cell);
		ReplaceComputationalCell(cmin, cell);
		// Find maximum.minimum neighbor values
		size_t nloop = neighbors.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			ComputationalCell const& cell_temp = *neighbors[i];
			if (!skip_key.empty() && *safe_retrieve(cell_temp.stickers.begin(), tracerstickernames.sticker_names.begin(),
				tracerstickernames.sticker_names.end(), skip_key))
				continue;
			if (tess.GetOriginalIndex(edge_list[i]->neighbors.first) == tess.GetOriginalIndex(edge_list[i]->neighbors.second))
				continue;
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			for (size_t j = 0; j < cell_temp.tracers.size(); ++j)
			{
				cmax.tracers[j] = std::max(cmax.tracers[j], cell_temp.tracers[j]);
				cmin.tracers[j] = std::min(cmin.tracers[j], cell_temp.tracers[j]);
			}
		}
		ReplaceComputationalCell(maxdiff, cmax);
		maxdiff -= cell;
		ReplaceComputationalCell(mindiff, cmin);
		mindiff -= cell;
		// limit the slope
		ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(*edge_list[0]), cm);
		ComputationalCell dphi = centroid_val - cell;
		vector<double> psi(4 + cell.tracers.size(), 1);
		const size_t nedges = edge_list.size();
		for (size_t i = 0; i < nedges; ++i)
		{
			if (tess.GetOriginalIndex(edge_list[i]->neighbors.first) == tess.GetOriginalIndex(edge_list[i]->neighbors.second))
				continue;
			if (i > 0)
			{
				ReplaceComputationalCell(centroid_val, cell);
				interp2(centroid_val, slope, CalcCentroid(*edge_list[i]), cm);
				ReplaceComputationalCell(dphi, centroid_val);
				dphi -= cell;
			}
			// density
			if (std::abs(dphi.density) > 0.1*std::max(std::abs(maxdiff.density), std::abs(mindiff.density)) || centroid_val.density*cell.density < 0)
			{
				if (dphi.density > 1e-9*cell.density)
					psi[0] = std::min(psi[0], maxdiff.density / dphi.density);
				else
					if (dphi.density<-1e-9*cell.density)
						psi[0] = std::min(psi[0], mindiff.density / dphi.density);
			}
			// pressure
			if (std::abs(dphi.pressure) > 0.1*std::max(std::abs(maxdiff.pressure), std::abs(mindiff.pressure)) || centroid_val.pressure*cell.pressure < 0)
			{
				if (dphi.pressure > 1e-9*cell.pressure)
					psi[1] = std::min(psi[1], maxdiff.pressure / dphi.pressure);
				else
					if (dphi.pressure<-1e-9*cell.pressure)
						psi[1] = std::min(psi[1], mindiff.pressure / dphi.pressure);
			}
			// xvelocity
			if (std::abs(dphi.velocity.x) > 0.1*std::max(std::abs(maxdiff.velocity.x), std::abs(mindiff.velocity.x)) || centroid_val.velocity.x*cell.velocity.x < 0)
			{
				if (dphi.velocity.x > std::abs(1e-9*cell.velocity.x))
					psi[2] = std::min(psi[2], maxdiff.velocity.x / dphi.velocity.x);
				else
					if (dphi.velocity.x<-std::abs(1e-9*cell.velocity.x))
						psi[2] = std::min(psi[2], mindiff.velocity.x / dphi.velocity.x);
			}
			// yvelocity
			if (std::abs(dphi.velocity.y) > 0.1*std::max(std::abs(maxdiff.velocity.y), std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (dphi.velocity.y > std::abs(1e-9*cell.velocity.y))
					psi[3] = std::min(psi[3], maxdiff.velocity.y / dphi.velocity.y);
				else
					if (dphi.velocity.y < -std::abs(1e-9*cell.velocity.y))
						psi[3] = std::min(psi[3], mindiff.velocity.y / dphi.velocity.y);
			}
			// tracers
			for (size_t j = 0; j < dphi.tracers.size(); ++j)
			{
				double cell_tracer = cell.tracers[j];
				double diff_tracer = maxdiff.tracers[j];
				if (std::abs(dphi.tracers[j]) > 0.1*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers[j])) || (
					centroid_val.tracers[j] * cell_tracer < 0))
				{
					if (dphi.tracers[j] > std::abs(1e-9*cell_tracer))
						psi[4 + j] = std::min(psi[4 + j], diff_tracer / dphi.tracers[j]);
					else
						if (dphi.tracers[j] < -std::abs(1e-9 * cell_tracer))
							psi[4 + j] = std::min(psi[4 + j], mindiff.tracers[j] / dphi.tracers[j]);
				}
			}
		}
		slope.xderivative.density *= psi[0];
		slope.yderivative.density *= psi[0];
		slope.xderivative.pressure *= psi[1];
		slope.yderivative.pressure *= psi[1];
		slope.xderivative.velocity.x *= psi[2];
		slope.yderivative.velocity.x *= psi[2];
		slope.xderivative.velocity.y *= psi[3];
		slope.yderivative.velocity.y *= psi[3];
		size_t counter = 4;
		size_t N = slope.xderivative.tracers.size();
		for (size_t k = 0; k < N; ++k)
		{
			slope.xderivative.tracers[k] *= psi[counter];
			slope.yderivative.tracers[k] *= psi[counter];
			++counter;
		}
	}

	void shocked_slope_limit(ComputationalCell const& cell, Vector2D const& cm,
		vector<ComputationalCell const*> const& neighbors, vector<const Edge*> const& edge_list,
		Slope  &slope, double diffusecoeff,TracerStickerNames const& tracerstickernames,
		string const& skip_key)
	{
		ComputationalCell cmax(cell), cmin(cell);
		// Find maximum values
		for (size_t i = 0; i < neighbors.size(); ++i)
		{
			ComputationalCell const& cell_temp = *neighbors[i];
			if (!skip_key.empty() && *safe_retrieve(cell_temp.stickers.begin(), tracerstickernames.sticker_names.begin(),
				tracerstickernames.sticker_names.end(),skip_key))
				continue;
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			for (size_t j = 0; j < cell_temp.tracers.size(); ++j)
			{
				cmax.tracers[j] = std::max(cmax.tracers[j], cell_temp.tracers[j]);
				cmin.tracers[j] = std::min(cmin.tracers[j], cell_temp.tracers[j]);
			}
		}
		ComputationalCell maxdiff = cmax - cell, mindiff = cmin - cell;
		// limit the slope
		vector<double> psi(4 + cell.tracers.size(), 1);
		for (size_t i = 0; i<edge_list.size(); ++i)
		{
			if (!skip_key.empty() && *safe_retrieve(neighbors[i]->stickers.begin(), 
				tracerstickernames.sticker_names.begin(), tracerstickernames.sticker_names.end(), skip_key))
				continue;
			ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(*edge_list[i]), cm);
			ComputationalCell dphi = centroid_val - cell;
			// density
			if (std::abs(dphi.density) > 0.1*std::max(std::abs(maxdiff.density), std::abs(mindiff.density)) || centroid_val.density*cell.density < 0)
			{
				if (std::abs(dphi.density) > 1e-9*cell.density)
					psi[0] = std::min(psi[0], std::max(diffusecoeff*(neighbors[i]->density - cell.density) / dphi.density, 0.0));
			}
			// pressure
			if (std::abs(dphi.pressure) > 0.1*std::max(std::abs(maxdiff.pressure), std::abs(mindiff.pressure)) || centroid_val.pressure*cell.pressure < 0)
			{
				if (std::abs(dphi.pressure) > 1e-9*cell.pressure)
					psi[1] = std::min(psi[1], std::max(diffusecoeff*(neighbors[i]->pressure - cell.pressure) / dphi.pressure, 0.0));
			}
			// xvelocity
			if (std::abs(dphi.velocity.x) > 0.1*std::max(std::abs(maxdiff.velocity.x), std::abs(mindiff.velocity.x)) || centroid_val.velocity.x*cell.velocity.x < 0)
			{
				if (std::abs(dphi.velocity.x) > 1e-9*cell.velocity.x)
					psi[2] = std::min(psi[2], std::max(diffusecoeff*(neighbors[i]->velocity.x - cell.velocity.x) / dphi.velocity.x, 0.0));
			}
			// yvelocity
			if (std::abs(dphi.velocity.y) > 0.1*std::max(std::abs(maxdiff.velocity.y), std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (std::abs(dphi.velocity.y) > 1e-9*cell.velocity.y)
					psi[3] = std::min(psi[3], std::max(diffusecoeff*(neighbors[i]->velocity.y - cell.velocity.y) / dphi.velocity.y, 0.0));
			}
			// tracers
			size_t counter = 0;
			for (size_t j = 0; j < dphi.tracers.size(); ++j)
			{
				double cell_tracer = cell.tracers[j];
				double diff_tracer = maxdiff.tracers[j];
				double centroid_tracer = centroid_val.tracers[j];
				if (std::abs(dphi.tracers[j]) > 0.1*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers[j])) ||
					centroid_tracer*cell_tracer < 0)
				{
					if (std::abs(dphi.tracers[j]) > std::abs(1e-9*cell_tracer))
						psi[4 + counter] = std::min(psi[4 + counter],
							std::max(diffusecoeff*(neighbors[i]->tracers[j] - cell_tracer) / dphi.tracers[j], 0.0));
				}
				++counter;
			}
		}
		slope.xderivative.density *= psi[0];
		slope.yderivative.density *= psi[0];
		slope.xderivative.pressure *= psi[1];
		slope.yderivative.pressure *= psi[1];
		slope.xderivative.velocity.x *= psi[2];
		slope.yderivative.velocity.x *= psi[2];
		slope.xderivative.velocity.y *= psi[3];
		slope.yderivative.velocity.y *= psi[3];
		size_t counter = 0;
		for (size_t k = 0; k < slope.xderivative.tracers.size(); ++k)
		{
			slope.xderivative.tracers[k] *= psi[4 + counter];
			slope.yderivative.tracers[k] *= psi[4 + counter];
			++counter;
		}
	}

	void GetBoundarySlope(ComputationalCell const& cell, Vector2D const& cell_cm,
		vector<ComputationalCell const*> const& neighbors, vector<Vector2D> const& neigh_cm,
		Slope &res)
	{
		size_t Nneigh = neigh_cm.size();
		ComputationalCell PhiSy, PhiSx;
	//	PhiSy.tracers.resize(cell.tracers.size(), 0);
	//	PhiSx.tracers.resize(cell.tracers.size(), 0);
		double SxSy(0), Sy2(0), Sx2(0);
		for (size_t i = 0; i < Nneigh; ++i)
		{
			PhiSy += (*neighbors[i] - cell)*(neigh_cm[i].y - cell_cm.y);
			PhiSx += (*neighbors[i] - cell)*(neigh_cm[i].x - cell_cm.x);
			SxSy += (neigh_cm[i].y - cell_cm.y)*(neigh_cm[i].x - cell_cm.x);
			Sx2 += (neigh_cm[i].x - cell_cm.x)*(neigh_cm[i].x - cell_cm.x);
			Sy2 += (neigh_cm[i].y - cell_cm.y)*(neigh_cm[i].y - cell_cm.y);
		}
		res.xderivative = (PhiSy*SxSy - PhiSx*Sy2) / (SxSy*SxSy - Sx2*Sy2);
		res.xderivative.stickers = cell.stickers;
		res.yderivative = (PhiSx*SxSy - PhiSy*Sx2) / (SxSy*SxSy - Sx2*Sy2);
		res.yderivative.stickers = cell.stickers;
	}


	void calc_slope
		(Tessellation const& tess,
			vector<ComputationalCell> const& cells,
			size_t cell_index,
			bool slf,
			double shockratio,
			double diffusecoeff,
			double pressure_ratio,
			EquationOfState const& eos,
			const vector<string>& flat_tracers,
			Slope &naive_slope_,
			Slope & res,
			Slope & temp1,
			ComputationalCell &temp2,
			ComputationalCell &temp3,
			ComputationalCell &temp4,
			ComputationalCell &temp5,
			vector<const Edge *> const& edge_list,
			vector<Vector2D> &neighbor_mesh_list,
			vector<Vector2D> &neighbor_cm_list,
			TracerStickerNames const& tracerstickernames,
			string const& skip_key)
	{
		GetNeighborMesh(tess, edge_list, cell_index, neighbor_mesh_list);
		GetNeighborCM(tess, edge_list, cell_index, neighbor_cm_list);
		vector<ComputationalCell const* > neighbor_list = GetNeighborCells(edge_list, cell_index, cells);

		ComputationalCell const& cell = cells[cell_index];
		bool boundary_slope = false;
		size_t Nneigh = neighbor_list.size();
		for (size_t i = 0; i < Nneigh; ++i)
			if (tess.GetOriginalIndex(edge_list[i]->neighbors.first) == tess.GetOriginalIndex(edge_list[i]->neighbors.second))
			{
				boundary_slope = true;
				break;
			}
		if (boundary_slope)
			GetBoundarySlope(cell, tess.GetCellCM(static_cast<int>(cell_index)),
				neighbor_list, neighbor_cm_list, res);
		else
			calc_naive_slope(cell, tess.GetMeshPoint(static_cast<int>(cell_index)), tess.GetCellCM(static_cast<int>(cell_index)),
				tess.GetVolume(static_cast<int>(cell_index)), neighbor_list, neighbor_mesh_list, neighbor_cm_list, edge_list,
				res, temp1);

		naive_slope_ = res;

		for (size_t i = 0; i < flat_tracers.size(); ++i)
		{
			size_t tindex = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
				flat_tracers[i]) - tracerstickernames.tracer_names.begin());
			assert(tindex < tracerstickernames.tracer_names.size());
			res.xderivative.tracers[tindex] = 0;
			res.yderivative.tracers[tindex] = 0;
		}

		if (slf)
		{
			if (!is_shock(res, tess.GetWidth(static_cast<int>(cell_index)), shockratio, cell, neighbor_list, pressure_ratio,
				eos.dp2c(cell.density, cell.pressure, cell.tracers,tracerstickernames.tracer_names)))
			{
				slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, res, temp2, temp3,
					temp4, temp5,tracerstickernames,skip_key,tess);
			}
			else
			{
				shocked_slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, res, diffusecoeff,tracerstickernames,skip_key);
			}
		}
	}
}

ComputationalCell LinearGaussImproved::Interp(ComputationalCell const& cell, size_t cell_index, Vector2D const& cm,
	Vector2D const& target)const
{
	return interp(cell, rslopes_[cell_index], target, cm);
}

LinearGaussImproved::LinearGaussImproved
(EquationOfState const& eos,
	GhostPointGenerator const& ghost,
	bool slf,
	double delta_v,
	double theta,
	double delta_P,
	const vector<string>& flat_tracers,
	string skip_key) :
	eos_(eos),
	ghost_(ghost),
	rslopes_(),
	naive_rslopes_(),
	slf_(slf),
	shockratio_(delta_v),
	diffusecoeff_(theta),
	pressure_ratio_(delta_P),
	flat_tracers_(flat_tracers),
	skip_key_(skip_key),
	to_skip_(){}

#ifdef RICH_MPI
namespace
{
	void exchange_ghost_slopes(Tessellation const& tess, vector<Slope> & slopes)
	{
		MPI_exchange_data(tess, slopes, true);
	}
}
#endif//RICH_MPI

void LinearGaussImproved::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells, double time, vector<pair<ComputationalCell, ComputationalCell> > &res,
	TracerStickerNames const& tracerstikersnames,CacheData const& /*cd*/) const
{
	const size_t CellNumber = static_cast<size_t>(tess.GetPointNo());
	vector<int> boundaryedges;
	// Get ghost points
	boost::container::flat_map<size_t, ComputationalCell> ghost_cells = ghost_.operator()(tess, cells, time, tracerstikersnames);
	// Copy ghost data into new cells vector
	vector<ComputationalCell> new_cells(cells);
	new_cells.resize
	  (static_cast<size_t>(tess.GetTotalPointNumber()));
	for (boost::container::flat_map<size_t, ComputationalCell>::const_iterator it = ghost_cells.begin(); it !=
		ghost_cells.end(); ++it)
		new_cells[it->first] = it->second;
	// Prepare slopes
	rslopes_.resize(CellNumber, Slope(cells[0], cells[0]));
	naive_rslopes_.resize(CellNumber);
	Slope temp1(cells[0], cells[0]);
	ComputationalCell temp2(cells[0]);
	ComputationalCell temp3(cells[0]);
	ComputationalCell temp4(cells[0]);
	ComputationalCell temp5(cells[0]);
	vector<const Edge *> edge_list;
	vector<Vector2D> neighbor_mesh_list;
	vector<Vector2D> neighbor_cm_list;
	for (size_t i = 0; i < CellNumber; ++i)
	{
		vector<int> const& edge_index = tess.GetCellEdges(static_cast<int>(i));
		GetEdgeList(tess, edge_index, edge_list);
		calc_slope(tess, new_cells, i, slf_, shockratio_, diffusecoeff_, pressure_ratio_, eos_,
			flat_tracers_, naive_rslopes_[i], rslopes_[i], temp1, temp2, temp3, temp4, temp5, edge_list,
			neighbor_mesh_list, neighbor_cm_list, tracerstikersnames,skip_key_);
		const size_t nloop = edge_index.size();
		for (size_t j = 0; j < nloop; ++j)
		{
			Edge const& edge = *edge_list[j];
			if (edge.neighbors.first == static_cast<int>(i))
			{
				ReplaceComputationalCell(res[static_cast<size_t>(edge_index[j])].first,
					new_cells[i]);
				interp2(res[static_cast<size_t>(edge_index[j])].first,
					rslopes_[i], CalcCentroid(edge), tess.GetCellCM(static_cast<int>(i)));
				if (edge.neighbors.second > static_cast<int>(CellNumber))
					boundaryedges.push_back(edge_index[j]);
			}
			else
			{
				ReplaceComputationalCell(res[static_cast<size_t>(edge_index[j])].second,
					new_cells[i]);
				interp2(res[static_cast<size_t>(edge_index[j])].second,
					rslopes_[i], CalcCentroid(edge), tess.GetCellCM(static_cast<int>(i)));
				if (edge.neighbors.first > static_cast<int>(CellNumber))
					boundaryedges.push_back(edge_index[j]);
			}
		}
	}
#ifdef RICH_MPI
	// communicate ghost slopes
	exchange_ghost_slopes(tess, rslopes_);
#endif //RICH_MPI
	// Interpolate the boundary edges
	for (size_t i = 0; i < boundaryedges.size(); ++i)
	{
		Edge const& edge = tess.GetEdge(boundaryedges[i]);
		if (edge.neighbors.first > static_cast<int>(CellNumber))
		{
			res[static_cast<size_t>(boundaryedges[i])].first = new_cells[static_cast<size_t>(edge.neighbors.first)];
#ifdef RICH_MPI
			if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
				interp2(res[static_cast<size_t>(boundaryedges[i])].first,
					rslopes_[static_cast<size_t>(edge.neighbors.first)], CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
			else
				res[static_cast<size_t>(boundaryedges[i])].first = interp(res[static_cast<size_t>(boundaryedges[i])].first,
					ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
						edge.neighbors.first), time, edge, tracerstikersnames), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
#else
			res[static_cast<size_t>(boundaryedges[i])].first = interp(res[static_cast<size_t>(boundaryedges[i])].first,
				ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
					edge.neighbors.first), time, edge, tracerstikersnames), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
#endif //RICH_MPI
		}
		else
		{
			res[static_cast<size_t>(boundaryedges[i])].second = new_cells[static_cast<size_t>(edge.neighbors.second)];
#ifdef RICH_MPI
			if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
				interp2(res[static_cast<size_t>(boundaryedges[i])].second,
					rslopes_[static_cast<size_t>(edge.neighbors.second)], CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
			else
				res[static_cast<size_t>(boundaryedges[i])].second = interp(res[static_cast<size_t>(boundaryedges[i])].second,
					ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
						edge.neighbors.second), time, edge, tracerstikersnames), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
#else
			res[static_cast<size_t>(boundaryedges[i])].second = interp(res[static_cast<size_t>(boundaryedges[i])].second,
				ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
					edge.neighbors.second), time, edge, tracerstikersnames), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
#endif //RICH_MPI
		}
	}
}


vector<Slope>& LinearGaussImproved::GetSlopes(void)const
{
	return rslopes_;
}

vector<Slope>& LinearGaussImproved::GetSlopesUnlimited(void)const
{
	return naive_rslopes_;
}

