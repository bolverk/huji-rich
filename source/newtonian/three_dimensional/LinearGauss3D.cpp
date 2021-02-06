#include "LinearGauss3D.hpp"
#include "../../misc/utils.hpp"
#include <array>
#include <iostream>
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

namespace
{
	void CheckCell(ComputationalCell3D const& cell)
	{
		if ((!(cell.density > 0)) || (!(cell.internal_energy > 0)) || (!std::isfinite(cell.velocity.x)) || (!std::isfinite(cell.velocity.y)) || (!std::isfinite(cell.velocity.z)))
			throw UniversalError("Bad cell after interpolation in LinearGauss3D");
	}

	void GetNeighborMesh(Tessellation3D const& tess, size_t cell_index,
		vector<Vector3D> &res, face_vec const& faces)
	{
		res.resize(faces.size());
		const int nloop = static_cast<int>(res.size());
		std::pair<size_t, size_t> neigh;
		for (int i = 0; i < nloop; ++i)
		{
			neigh = tess.GetFaceNeighbors(faces[i]);
			if (neigh.first == cell_index)
				res[i] = tess.GetMeshPoint(neigh.second);
			else
				res[i] = tess.GetMeshPoint(neigh.first);
		}
	}

	void GetNeighborCM(Tessellation3D const& tess, size_t cell_index,
		vector<Vector3D> &res, face_vec const& faces)
	{
		res.resize(faces.size());
		const size_t nloop = faces.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			if (tess.GetFaceNeighbors(faces[i]).first == cell_index)
				res[i] = tess.GetCellCM(tess.GetFaceNeighbors(faces[i]).second);
			else
				res[i] = tess.GetCellCM(tess.GetFaceNeighbors(faces[i]).first);
		}
	}

	void GetNeighborCells(Tessellation3D const& tess, size_t cell_index,
		vector<ComputationalCell3D> const& cells, face_vec const& faces, vector<ComputationalCell3D> &res)
	{
		const size_t nloop = faces.size();
		res.resize(nloop);
		for (size_t i = 0; i < nloop; ++i)
		{
			size_t other_cell = (tess.GetFaceNeighbors(faces[i]).first == cell_index) ?
				tess.GetFaceNeighbors(faces[i]).second : tess.GetFaceNeighbors(faces[i]).first;
			ReplaceComputationalCell(res[i], cells[other_cell]);
		}
	}

	void calc_naive_slope(ComputationalCell3D const& cell,
		Vector3D const& center, Vector3D const& cell_cm, double cell_volume, vector<ComputationalCell3D> const& neighbors,
		vector<Vector3D> const& neighbor_centers,
		vector<Vector3D> const& neigh_cm, Tessellation3D const& tess,
		Slope3D &res, Slope3D &temp, size_t /*index*/, face_vec const& faces,
		std::vector<Vector3D> c_ij)
	{
		size_t n = neighbor_centers.size();
		if (n > 60)
			std::cout << "Cell has too many neighbors in calc naive slope, Cell x cor " << center.x <<
			" Cell y cor " << center.y << " Cell z cor " << center.z << std::endl;
		// Create the matrix to invert and the vector to compare
		std::array<double, 9>  m;
		std::fill_n(m.begin(), 9, 0.0);
		c_ij.resize(n);

		for (size_t i = 0; i < n; i++)
		{
			c_ij[i] = neigh_cm[i];
			c_ij[i] += cell_cm;
			c_ij[i] *= -0.5;
			c_ij[i] += tess.FaceCM(faces[i]);
			const Vector3D r_ij = normalize(neighbor_centers[i] - center);
			const double A = tess.GetArea(faces[i]);
			m[0] -= c_ij[i].x*r_ij.x*A;
			m[1] -= c_ij[i].y*r_ij.x*A;
			m[2] -= c_ij[i].z*r_ij.x*A;
			m[3] -= c_ij[i].x*r_ij.y*A;
			m[4] -= c_ij[i].y*r_ij.y*A;
			m[5] -= c_ij[i].z*r_ij.y*A;
			m[6] -= c_ij[i].x*r_ij.z*A;
			m[7] -= c_ij[i].y*r_ij.z*A;
			m[8] -= c_ij[i].z*r_ij.z*A;
			if (i == 0)
			{
				ReplaceComputationalCell(temp.xderivative, neighbors[i]);
				temp.xderivative *= r_ij.x*A;
				ReplaceComputationalCell(temp.yderivative, neighbors[i]);
				temp.yderivative *= r_ij.y*A;
				ReplaceComputationalCell(temp.zderivative, neighbors[i]);
				temp.zderivative *= r_ij.z*A;
			}
			else
			{
				ComputationalCellAddMult(temp.xderivative, neighbors[i], r_ij.x*A);
				ComputationalCellAddMult(temp.yderivative, neighbors[i], r_ij.y*A);
				ComputationalCellAddMult(temp.zderivative, neighbors[i], r_ij.z*A);
			}
			ComputationalCellAddMult(temp.xderivative, cell, r_ij.x*A);
			ComputationalCellAddMult(temp.yderivative, cell, r_ij.y*A);
			ComputationalCellAddMult(temp.zderivative, cell, r_ij.z*A);
		}
		double v_inv = 1.0 / cell_volume;
		for (size_t i = 0; i < 9; ++i)
			m[i] *= v_inv;
		m[0] += 1;
		m[4] += 1;
		m[8] += 1;
		// Find the det
		const double det = -m[2] * m[4] * m[6] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[0] * m[5] * m[7] -
			m[1] * m[3] * m[8] + m[0] * m[4] * m[8];
		// Check none singular
		if (std::abs(det) < 1e-10)
		{
			UniversalError eo("Singular matrix");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			eo.AddEntry("Cell z cor", center.z);
			eo.AddEntry("Cell CMx cor", cell_cm.x);
			eo.AddEntry("Cell CMy cor", cell_cm.y);
			eo.AddEntry("Cell CMz cor", cell_cm.z);
			eo.AddEntry("Cell volume", cell_volume);
			eo.AddEntry("Det was", det);
			for (size_t i = 0; i < faces.size(); ++i)
			{
				c_ij[0] = tess.FaceCM(faces[i]) - 0.5 * (neigh_cm[i] + cell_cm);
				eo.AddEntry("Neighbor x", neighbor_centers[i].x);
				eo.AddEntry("Neighbor y", neighbor_centers[i].y);
				eo.AddEntry("Neighbor z", neighbor_centers[i].z);
				eo.AddEntry("Face", static_cast<double>(faces[i]));
				eo.AddEntry("Neighbor Cx", c_ij[0].x);
				eo.AddEntry("Neighbor Cy", c_ij[0].y);
				eo.AddEntry("Neighbor Cz", c_ij[0].z);
				eo.AddEntry("Face Cx", tess.FaceCM(faces[i]).x);
				eo.AddEntry("Face Cy", tess.FaceCM(faces[i]).y);
				eo.AddEntry("Face Cz", tess.FaceCM(faces[i]).z);
				eo.AddEntry("Face area", tess.GetArea(faces[i]));
			}
			for (size_t i = 0; i < 9; ++i)
				eo.AddEntry("M", m[i]);
			throw eo;
		}
		// Invert the matrix
		std::array<double, 9>  m_inv;
		std::fill_n(m_inv.begin(), 9, 0);
		m_inv[0] = m[4] * m[8] - m[5] * m[7];
		m_inv[1] = m[2] * m[7] - m[1] * m[8];
		m_inv[2] = m[1] * m[5] - m[2] * m[4];
		m_inv[3] = m[5] * m[6] - m[3] * m[8];
		m_inv[4] = m[0] * m[8] - m[2] * m[6];
		m_inv[5] = m[2] * m[3] - m[5] * m[0];
		m_inv[6] = m[3] * m[7] - m[6] * m[4];
		m_inv[7] = m[6] * m[1] - m[0] * m[7];
		m_inv[8] = m[4] * m[0] - m[1] * m[3];
		for (size_t i = 0; i < 9; ++i)
			m_inv[i] /= (2 * cell_volume * det);

		ReplaceComputationalCell(res.xderivative, temp.xderivative);
		res.xderivative *= m_inv[0];
		ComputationalCellAddMult(res.xderivative, temp.yderivative, m_inv[1]);
		ComputationalCellAddMult(res.xderivative, temp.zderivative, m_inv[2]);

		ReplaceComputationalCell(res.yderivative, temp.xderivative);
		res.yderivative *= m_inv[3];
		ComputationalCellAddMult(res.yderivative, temp.yderivative, m_inv[4]);
		ComputationalCellAddMult(res.yderivative, temp.zderivative, m_inv[5]);

		ReplaceComputationalCell(res.zderivative, temp.xderivative);
		res.zderivative *= m_inv[6];
		ComputationalCellAddMult(res.zderivative, temp.yderivative, m_inv[7]);
		ComputationalCellAddMult(res.zderivative, temp.zderivative, m_inv[8]);
	}


	double PressureRatio(ComputationalCell3D const& cell, vector<ComputationalCell3D> const& neigh)
	{
		double res = 1.0;
		double p = cell.pressure;
		size_t N = neigh.size();
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(min:res)
#endif
		for (size_t i = 0; i < N; i++)
		{
			if (p > neigh[i].pressure)
				res = std::min(res, neigh[i].pressure / p);
			else
				res = std::min(res, p / neigh[i].pressure);
		}
		return res;
	}

	bool is_shock(Slope3D const& naive_slope, double cell_width, double shock_ratio,
		ComputationalCell3D const& cell, vector<ComputationalCell3D> const& neighbor_list, double pressure_ratio, double cs)
	{
		const bool cond1 = (naive_slope.xderivative.velocity.x + naive_slope.yderivative.velocity.y + naive_slope.zderivative.velocity.z)*
			cell_width < (-shock_ratio)*cs;
		const bool cond2 = PressureRatio(cell, neighbor_list) < pressure_ratio;
		return cond1 || cond2;
	}

	ComputationalCell3D interp(ComputationalCell3D const& cell, Slope3D const& slope,
		Vector3D const& target, Vector3D const& cm, EquationOfState const& eos, TracerStickerNames const& tsn,
		bool pressure_calc)
	{
		ComputationalCell3D res(cell);
		ComputationalCellAddMult(res, slope.xderivative, target.x - cm.x);
		ComputationalCellAddMult(res, slope.yderivative, target.y - cm.y);
		ComputationalCellAddMult(res, slope.zderivative, target.z - cm.z);
		if (pressure_calc)
			try
		{
			//res.pressure = eos.de2p(res.density, res.internal_energy, res.tracers, tsn.tracer_names);
			res.internal_energy = eos.dp2e(res.density, res.pressure, res.tracers, tsn.tracer_names);
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("density", res.density);
			eo.AddEntry("internal energy", res.internal_energy);
			throw eo;
		}
		return res;
	}

	void interp23Dsimple(ComputationalCell3D &res, Slope3D const& slope,
		Vector3D const& target, Vector3D const& cm)
	{
		ComputationalCellAddMult(res, slope.xderivative, target.x - cm.x);
		ComputationalCellAddMult(res, slope.yderivative, target.y - cm.y);
		ComputationalCellAddMult(res, slope.zderivative, target.z - cm.z);
	}

	void interp23D(ComputationalCell3D &res, Slope3D const& slope,
		Vector3D const& target, Vector3D const& cm, EquationOfState const& eos, TracerStickerNames const& tsn,
		bool pressure_calc)
	{
		ComputationalCellAddMult(res, slope.xderivative, target.x - cm.x);
		ComputationalCellAddMult(res, slope.yderivative, target.y - cm.y);
		ComputationalCellAddMult(res, slope.zderivative, target.z - cm.z);
		try
		{
			if (!pressure_calc)
				res.pressure = eos.de2p(res.density, res.internal_energy, res.tracers, tsn.tracer_names);
			else
				res.internal_energy = eos.dp2e(res.density, res.pressure, res.tracers, tsn.tracer_names);
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("density", res.density);
			eo.AddEntry("internal energy", res.internal_energy);
			throw eo;
		}
	}

	void slope_limit(ComputationalCell3D const& cell, Vector3D const& cm,
		vector<ComputationalCell3D> const& neighbors, Slope3D &slope, ComputationalCell3D &cmax,
		ComputationalCell3D &cmin, ComputationalCell3D &maxdiff, ComputationalCell3D &mindiff,
		TracerStickerNames const& tracerstickernames, string const& skip_key, Tessellation3D const& tess,
		size_t /*cell_index*/, face_vec const& faces, EquationOfState const& eos)
	{
		ReplaceComputationalCell(cmax, cell);
		ReplaceComputationalCell(cmin, cell);
		// Find maximum.minimum neighbor values
		size_t nloop = neighbors.size();
		size_t ntracer = cell.tracers.size();
		for (size_t i = 0; i < nloop; ++i)
		{
			ComputationalCell3D const& cell_temp = neighbors[i];
			if (!skip_key.empty() && *safe_retrieve(cell_temp.stickers.begin(), tracerstickernames.sticker_names.begin(),
				tracerstickernames.sticker_names.end(), skip_key))
				continue;
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmax.velocity.z = std::max(cmax.velocity.z, cell_temp.velocity.z);
			cmax.internal_energy = std::max(cmax.internal_energy, cell_temp.internal_energy);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			cmin.velocity.z = std::min(cmin.velocity.z, cell_temp.velocity.z);
			cmin.internal_energy = std::min(cmin.internal_energy, cell_temp.internal_energy);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
			for (size_t j = 0; j < ntracer; ++j)
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
		ComputationalCell3D centroid_val = interp(cell, slope, tess.FaceCM(faces[0]), cm, eos, tracerstickernames, false);
		ComputationalCell3D dphi = centroid_val - cell;
		vector<double> psi(6 + cell.tracers.size(), 1);
		const size_t nedges = faces.size();
		const double skipfactor = 1e-3;
		for (size_t i = 0; i < nedges; i++)
		{
			if (i > 0)
			{
				ReplaceComputationalCell(centroid_val, cell);
				interp23Dsimple(centroid_val, slope, tess.FaceCM(faces[i]), cm);
				ReplaceComputationalCell(dphi, centroid_val);
				dphi -= cell;
			}
			// density
			if (std::abs(dphi.density) > skipfactor*std::max(std::abs(maxdiff.density), std::abs(mindiff.density)) || centroid_val.density*cell.density < 0)
			{
				if (dphi.density > 1e-9*cell.density)
					psi[0] = std::min(psi[0], std::max(maxdiff.density / dphi.density, 0.0));
				else
					if (dphi.density < -1e-9*cell.density)
						psi[0] = std::min(psi[0], std::max(mindiff.density / dphi.density, 0.0));
			}
			// pressure
			if (std::abs(dphi.pressure) > skipfactor*std::max(std::abs(maxdiff.pressure), std::abs(mindiff.pressure)) || centroid_val.pressure*cell.pressure < 0)
			{
				if (dphi.pressure > 1e-9*cell.pressure)
					psi[1] = std::min(psi[1], std::max(maxdiff.pressure / dphi.pressure, 0.0));
				else
					if (dphi.pressure < -1e-9*cell.pressure)
						psi[1] = std::min(psi[1], std::max(mindiff.pressure / dphi.pressure, 0.0));
			}
			// internal_energy
			if (std::abs(dphi.internal_energy) > skipfactor*std::max(std::abs(maxdiff.internal_energy), std::abs(mindiff.internal_energy)) || centroid_val.internal_energy*cell.internal_energy < 0)
			{
				if (dphi.internal_energy > 1e-9*cell.internal_energy)
					psi[5] = std::min(psi[5], std::max(maxdiff.internal_energy / dphi.internal_energy, 0.0));
				else
					if (dphi.internal_energy < -1e-9*cell.internal_energy)
						psi[5] = std::min(psi[5], std::max(mindiff.internal_energy / dphi.internal_energy, 0.0));
			}
			// xvelocity
			if (std::abs(dphi.velocity.x) > skipfactor*std::max(std::abs(maxdiff.velocity.x), std::abs(mindiff.velocity.x)) || centroid_val.velocity.x*cell.velocity.x < 0)
			{
				if (dphi.velocity.x > std::abs(1e-9*cell.velocity.x))
					psi[2] = std::min(psi[2], std::max(maxdiff.velocity.x / dphi.velocity.x, 0.0));
				else
					if (dphi.velocity.x < -std::abs(1e-9*cell.velocity.x))
						psi[2] = std::min(psi[2], std::max(mindiff.velocity.x / dphi.velocity.x, 0.0));
			}
			// yvelocity
			if (std::abs(dphi.velocity.y) > skipfactor*std::max(std::abs(maxdiff.velocity.y), std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (dphi.velocity.y > std::abs(1e-9*cell.velocity.y))
					psi[3] = std::min(psi[3], std::max(maxdiff.velocity.y / dphi.velocity.y, 0.0));
				else
					if (dphi.velocity.y < -std::abs(1e-9*cell.velocity.y))
						psi[3] = std::min(psi[3], std::max(mindiff.velocity.y / dphi.velocity.y, 0.0));
			}
			// zvelocity
			if (std::abs(dphi.velocity.z) > skipfactor*std::max(std::abs(maxdiff.velocity.z), std::abs(mindiff.velocity.z)) || centroid_val.velocity.z*cell.velocity.z < 0)
			{
				if (dphi.velocity.z > std::abs(1e-9*cell.velocity.z))
					psi[4] = std::min(psi[4], std::max(maxdiff.velocity.z / dphi.velocity.z, 0.0));
				else
					if (dphi.velocity.z < -std::abs(1e-9*cell.velocity.z))
						psi[4] = std::min(psi[4], std::max(mindiff.velocity.z / dphi.velocity.z, 0.0));
			}
			// tracers
			for (size_t j = 0; j < ntracer; ++j)
			{
				double cell_tracer = cell.tracers[j];
				double diff_tracer = maxdiff.tracers[j];
				if (std::abs(dphi.tracers[j]) > skipfactor*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers[j])) || (
					centroid_val.tracers[j] * cell_tracer < 0))
				{
					if (dphi.tracers[j] > std::abs(1e-9*cell_tracer))
						psi[6 + j] = std::min(psi[6 + j], std::max(diff_tracer / dphi.tracers[j], 0.0));
					else
						if (dphi.tracers[j] < -std::abs(1e-9 * cell_tracer))
							psi[6 + j] = std::min(psi[6 + j], std::max(mindiff.tracers[j] / dphi.tracers[j], 0.0));
				}
			}
		}
		psi[1] = std::min(psi[1], psi[5]);
		psi[5] = psi[1];

		slope.xderivative.density *= psi[0];
		slope.yderivative.density *= psi[0];
		slope.zderivative.density *= psi[0];
		slope.xderivative.pressure *= psi[1];
		slope.yderivative.pressure *= psi[1];
		slope.zderivative.pressure *= psi[1];
		slope.xderivative.velocity.x *= psi[2];
		slope.yderivative.velocity.x *= psi[2];
		slope.zderivative.velocity.x *= psi[2];
		slope.xderivative.velocity.y *= psi[3];
		slope.yderivative.velocity.y *= psi[3];
		slope.zderivative.velocity.y *= psi[3];
		slope.xderivative.velocity.z *= psi[4];
		slope.yderivative.velocity.z *= psi[4];
		slope.zderivative.velocity.z *= psi[4];
		slope.xderivative.internal_energy *= psi[5];
		slope.yderivative.internal_energy *= psi[5];
		slope.zderivative.internal_energy *= psi[5];
		size_t counter = 6;
		size_t N = slope.xderivative.tracers.size();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t k = 0; k < N; ++k)
		{
			slope.xderivative.tracers[k] *= psi[counter];
			slope.yderivative.tracers[k] *= psi[counter];
			slope.zderivative.tracers[k] *= psi[counter];
			++counter;
		}
	}

	void shocked_slope_limit(ComputationalCell3D const& cell, Vector3D const& cm,
		vector<ComputationalCell3D> const& neighbors,
		Slope3D  &slope, double diffusecoeff, TracerStickerNames const& tracerstickernames,
		string const& skip_key, Tessellation3D const& tess, size_t cell_index, face_vec const& faces,
		EquationOfState const& eos)
	{
		const double small_factor = 1e-9;
		ComputationalCell3D cmax(cell), cmin(cell);
		size_t N = faces.size();
		size_t ntracer = cell.tracers.size();
		// Find maximum values
		for (size_t i = 0; i < N; ++i)
		{
			ComputationalCell3D const& cell_temp = neighbors[i];
			if (!skip_key.empty() && *safe_retrieve(cell_temp.stickers.begin(), tracerstickernames.sticker_names.begin(),
				tracerstickernames.sticker_names.end(), skip_key))
				continue;
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.internal_energy = std::max(cmax.internal_energy, cell_temp.internal_energy);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmax.velocity.z = std::max(cmax.velocity.z, cell_temp.velocity.z);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.internal_energy = std::min(cmin.internal_energy, cell_temp.internal_energy);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			cmin.velocity.z = std::min(cmin.velocity.z, cell_temp.velocity.z);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
			for (size_t j = 0; j < ntracer; ++j)
			{
				cmax.tracers[j] = std::max(cmax.tracers[j], cell_temp.tracers[j]);
				cmin.tracers[j] = std::min(cmin.tracers[j], cell_temp.tracers[j]);
			}
		}
		ComputationalCell3D maxdiff = cmax - cell, mindiff = cmin - cell;
		// limit the slope
		vector<double> psi(6 + cell.tracers.size(), 1);
		for (size_t i = 0; i < N; ++i)
		{
			if (!skip_key.empty() && *safe_retrieve(neighbors[i].stickers.begin(),
				tracerstickernames.sticker_names.begin(), tracerstickernames.sticker_names.end(), skip_key))
				continue;
			ComputationalCell3D centroid_val = interp(cell, slope, tess.FaceCM(faces[i]), cm, eos, tracerstickernames, false);
			ComputationalCell3D dphi = centroid_val - cell;
			// density
			if (std::abs(dphi.density) > small_factor*std::max(std::abs(maxdiff.density), std::abs(mindiff.density)) || centroid_val.density*cell.density < 0)
			{
				if (std::abs(dphi.density) > 1e-9*cell.density)
					psi[0] = std::min(psi[0], std::max(diffusecoeff*(neighbors[i].density - cell.density) / dphi.density, 0.0));
			}
			// pressure
			if (std::abs(dphi.pressure) > small_factor*std::max(std::abs(maxdiff.pressure), std::abs(mindiff.pressure)) || centroid_val.pressure*cell.pressure < 0)
			{
				if (std::abs(dphi.pressure) > 1e-9*cell.pressure)
					psi[1] = std::min(psi[1], std::max(diffusecoeff*(neighbors[i].pressure - cell.pressure) / dphi.pressure, 0.0));
			}
			// internal_energy
			if (std::abs(dphi.internal_energy) > small_factor*std::max(std::abs(maxdiff.internal_energy), std::abs(mindiff.internal_energy)) || centroid_val.internal_energy*cell.internal_energy < 0)
			{
				if (std::abs(dphi.internal_energy) > 1e-9*cell.internal_energy)
					psi[5] = std::min(psi[5], std::max(diffusecoeff*(neighbors[i].internal_energy - cell.internal_energy) / dphi.internal_energy, 0.0));
			}
			// xvelocity
			if (std::abs(dphi.velocity.x) > small_factor*std::max(std::abs(maxdiff.velocity.x), std::abs(mindiff.velocity.x)) || centroid_val.velocity.x*cell.velocity.x < 0)
			{
				if (std::abs(dphi.velocity.x) > 1e-9*std::abs(cell.velocity.x))
					psi[2] = std::min(psi[2], std::max(diffusecoeff*(neighbors[i].velocity.x - cell.velocity.x) / dphi.velocity.x, 0.0));
			}
			// yvelocity
			if (std::abs(dphi.velocity.y) > small_factor*std::max(std::abs(maxdiff.velocity.y), std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (std::abs(dphi.velocity.y) > 1e-9*std::abs(cell.velocity.y))
					psi[3] = std::min(psi[3], std::max(diffusecoeff*(neighbors[i].velocity.y - cell.velocity.y) / dphi.velocity.y, 0.0));
			}
			// zvelocity
			if (std::abs(dphi.velocity.z) > small_factor*std::max(std::abs(maxdiff.velocity.z), std::abs(mindiff.velocity.z)) || centroid_val.velocity.z*cell.velocity.z < 0)
			{
				if (std::abs(dphi.velocity.z) > 1e-9*std::abs(cell.velocity.z))
					psi[4] = std::min(psi[4], std::max(diffusecoeff*(neighbors[i].velocity.z - cell.velocity.z) / dphi.velocity.z, 0.0));
			}
			// tracers
			size_t counter = 0;
			for (size_t j = 0; j < ntracer; ++j)
			{
				double cell_tracer = cell.tracers[j];
				double diff_tracer = maxdiff.tracers[j];
				double centroid_tracer = centroid_val.tracers[j];
				if (std::abs(dphi.tracers[j]) > 0.001*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers[j])) ||
					centroid_tracer * cell_tracer < 0)
				{
					if (std::abs(dphi.tracers[j]) > std::abs(1e-9*cell_tracer))
						psi[6 + counter] = std::min(psi[6 + counter],
							std::max(diffusecoeff*(neighbors[i].tracers[j] - cell_tracer) / dphi.tracers[j], 0.0));
				}
				++counter;
			}
		}

		psi[1] = std::min(psi[1], psi[5]);
		psi[5] = psi[1];

		slope.xderivative.density *= psi[0];
		slope.yderivative.density *= psi[0];
		slope.zderivative.density *= psi[0];
		slope.xderivative.pressure *= psi[1];
		slope.yderivative.pressure *= psi[1];
		slope.zderivative.pressure *= psi[1];
		slope.xderivative.velocity.x *= psi[2];
		slope.yderivative.velocity.x *= psi[2];
		slope.zderivative.velocity.x *= psi[2];
		slope.xderivative.velocity.y *= psi[3];
		slope.yderivative.velocity.y *= psi[3];
		slope.zderivative.velocity.y *= psi[3];
		slope.xderivative.velocity.z *= psi[4];
		slope.yderivative.velocity.z *= psi[4];
		slope.zderivative.velocity.z *= psi[4];
		slope.xderivative.internal_energy *= psi[5];
		slope.yderivative.internal_energy *= psi[5];
		slope.zderivative.internal_energy *= psi[5];
		size_t counter = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
		for (size_t k = 0; k < ntracer; ++k)
		{
			slope.xderivative.tracers[k] *= psi[6 + counter];
			slope.yderivative.tracers[k] *= psi[6 + counter];
			slope.zderivative.tracers[k] *= psi[6 + counter];
			++counter;
		}
		// make sure velocity slope is not too large in supersonic regions
		double maxDv = ScalarProd(slope.xderivative.velocity, slope.xderivative.velocity)
			+ ScalarProd(slope.yderivative.velocity, slope.yderivative.velocity) +
			ScalarProd(slope.zderivative.velocity, slope.zderivative.velocity);
		maxDv *= tess.GetWidth(cell_index) * tess.GetWidth(cell_index);
		if (maxDv > 100 * ScalarProd(cell.velocity, cell.velocity))
		{
			double sfactor = fastsqrt(100 * ScalarProd(cell.velocity, cell.velocity) / maxDv);
			slope.xderivative.velocity.x *= sfactor;
			slope.yderivative.velocity.x *= sfactor;
			slope.zderivative.velocity.x *= sfactor;
			slope.xderivative.velocity.y *= sfactor;
			slope.yderivative.velocity.y *= sfactor;
			slope.zderivative.velocity.y *= sfactor;
			slope.xderivative.velocity.z *= sfactor;
			slope.yderivative.velocity.z *= sfactor;
			slope.zderivative.velocity.z *= sfactor;
		}
	}

	void GetBoundarySlope(ComputationalCell3D const& cell, Vector3D const& cell_cm,
		vector<ComputationalCell3D> const& neighbors,
		vector<Vector3D> const& neigh_cm,
		Slope3D &res)
	{
		size_t Nneigh = neigh_cm.size();
		ComputationalCell3D PhiSy, PhiSx, PhiSz;
		//		PhiSy.tracers.resize(cell.tracers.size(), 0);
			//	PhiSx.tracers.resize(cell.tracers.size(), 0);
				//PhiSz.tracers.resize(cell.tracers.size(), 0);
		double SxSy(0), Sy2(0), Sx2(0), SxSz(0), Sz2(0), SzSy(0);
		for (size_t i = 0; i < Nneigh; ++i)
		{
			PhiSy += (neighbors[i] - cell)*(neigh_cm[i].y - cell_cm.y);
			PhiSx += (neighbors[i] - cell)*(neigh_cm[i].x - cell_cm.x);
			PhiSz += (neighbors[i] - cell)*(neigh_cm[i].z - cell_cm.z);
			SxSy += (neigh_cm[i].y - cell_cm.y)*(neigh_cm[i].x - cell_cm.x);
			Sx2 += (neigh_cm[i].x - cell_cm.x)*(neigh_cm[i].x - cell_cm.x);
			Sy2 += (neigh_cm[i].y - cell_cm.y)*(neigh_cm[i].y - cell_cm.y);
			SxSz += (neigh_cm[i].z - cell_cm.z)*(neigh_cm[i].x - cell_cm.x);
			SzSy += (neigh_cm[i].z - cell_cm.z)*(neigh_cm[i].y - cell_cm.y);
			Sz2 += (neigh_cm[i].z - cell_cm.z)*(neigh_cm[i].z - cell_cm.z);
		}
		double bottom = 1.0 / SxSz * SxSz*Sy2 + SxSy * SxSy*Sz2 - Sx2 * Sy2*Sz2 - 2 * SxSy*SxSz*SzSy + Sx2 * SzSy*Sz2;
		res.xderivative = (PhiSz*SxSz*Sy2 + PhiSy * SxSy*Sz2 - PhiSx * Sy2*Sz2 - PhiSz * SxSy*SzSy - PhiSy * SxSz*SzSy +
			PhiSx * SzSy*SzSy)*bottom;
		res.yderivative = (PhiSz*SzSy*Sx2 + PhiSy * SxSz*SxSz - PhiSx * SxSz*SzSy - PhiSz * SxSy*SxSz - PhiSy * Sx2*Sz2 +
			PhiSx * SxSy*Sz2)*bottom;
		res.zderivative = (PhiSz*SxSy*SxSy - PhiSy * SxSy*SxSz - PhiSz * Sy2*Sx2 + PhiSx * SxSz*Sy2 + PhiSy * Sx2*SzSy -
			PhiSx * SxSy*SzSy) *bottom;
		res.xderivative.stickers = cell.stickers;
		res.yderivative.stickers = cell.stickers;
		res.zderivative.stickers = cell.stickers;
	}


	void calc_slope(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, size_t cell_index, bool slf,
		double shockratio, double diffusecoeff, double pressure_ratio, EquationOfState const& eos,
		const vector<string>& calc_tracers, Slope3D &naive_slope_, Slope3D & res, Slope3D &temp1, ComputationalCell3D &temp2,
		ComputationalCell3D &temp3, ComputationalCell3D &temp4, ComputationalCell3D &temp5,
		vector<Vector3D> &neighbor_mesh_list,
		vector<Vector3D> &neighbor_cm_list,
		TracerStickerNames const& tracerstickernames, string const& skip_key,
		std::vector<Vector3D> &c_ij, vector<ComputationalCell3D> &neighbor_list)
	{
		face_vec const& faces = tess.GetCellFaces(cell_index);
		GetNeighborMesh(tess, cell_index, neighbor_mesh_list, faces);
		GetNeighborCM(tess, cell_index, neighbor_cm_list, faces);
		GetNeighborCells(tess, cell_index, cells, faces, neighbor_list);
		ComputationalCell3D const& cell = cells[cell_index];
		bool boundary_slope = false;
		size_t Nneigh = faces.size();
		for (size_t i = 0; i < Nneigh; ++i)
			if (tess.BoundaryFace(faces[i]))
			{
				boundary_slope = true;
				break;
			}
		if (boundary_slope)
			GetBoundarySlope(cell, tess.GetCellCM(cell_index), neighbor_list, neighbor_cm_list, res);
		else
			calc_naive_slope(cell, tess.GetMeshPoint(cell_index), tess.GetCellCM(cell_index),
				tess.GetVolume(cell_index), neighbor_list, neighbor_mesh_list, neighbor_cm_list, tess, res, temp1,
				cell_index, faces, c_ij);

		naive_slope_ = res;

		for (size_t i = 0; i < tracerstickernames.tracer_names.size(); ++i)
		{
			if (std::find(calc_tracers.begin(), calc_tracers.end(), tracerstickernames.tracer_names[i]) == calc_tracers.end())
			{
				res.xderivative.tracers[i] = 0;
				res.yderivative.tracers[i] = 0;
				res.zderivative.tracers[i] = 0;
			}
		}

		if (slf)
		{
#ifdef RICH_DEBUG
			try
			{
#endif
				if (!is_shock(res, tess.GetWidth(cell_index), shockratio, cell, neighbor_list, pressure_ratio,
					eos.dp2c(cell.density, cell.pressure, cell.tracers, tracerstickernames.tracer_names)))
				{
					slope_limit(cell, tess.GetCellCM(cell_index), neighbor_list, res, temp2, temp3, temp4, temp5,
						tracerstickernames, skip_key, tess, cell_index, faces, eos);
				}
				else
				{
					shocked_slope_limit(cell, tess.GetCellCM(cell_index), neighbor_list, res, diffusecoeff, tracerstickernames,
						skip_key, tess, cell_index, faces, eos);
				}
#ifdef RICH_DEBUG
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Error LinearGauss3D", 0);
				eo.AddEntry("Cell number", cell_index);
				throw eo;
			}
#endif
		}
	}

#ifdef RICH_MPI

	void exchange_ghost_slopes(Tessellation3D const& tess, vector<Slope3D> & slopes)
	{
		Slope3D sdummy;
		MPI_exchange_data(tess, slopes, true,&sdummy);
	}
#endif//RICH_MPI
}

void LinearGauss3D::Interp(ComputationalCell3D &res, ComputationalCell3D const& cell, size_t cell_index, Vector3D const& cm,
	Vector3D const& target, EquationOfState const& eos, TracerStickerNames const& tsn)const
{
	try
	{
		res = interp(cell, rslopes_[cell_index], target, cm, eos, tsn, true);
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Cell density", cell.density);
		eo.AddEntry("Cell internal energy", cell.internal_energy);
		eo.AddEntry("cell index", static_cast<double>(cell_index));
		eo.AddEntry("CMx", cm.x);
		eo.AddEntry("CMy", cm.y);
		eo.AddEntry("CMz", cm.z);
		eo.AddEntry("Targetx", target.x);
		eo.AddEntry("Targety", target.y);
		eo.AddEntry("Targetz", target.z);
		throw eo;
	}
}

LinearGauss3D::LinearGauss3D(EquationOfState const& eos, TracerStickerNames const& tsn, Ghost3D const& ghost, bool slf, double delta_v, double theta,
	double delta_P, bool SR, const vector<string>& calc_tracers, string skip_key,bool pressure_calc) : eos_(eos), tsn_(tsn), ghost_(ghost), rslopes_(),
	naive_rslopes_(), slf_(slf), shockratio_(delta_v), diffusecoeff_(theta), pressure_ratio_(delta_P), SR_(SR),
	calc_tracers_(calc_tracers), skip_key_(skip_key), to_skip_(),pressure_calc_(pressure_calc) {}

void LinearGauss3D::BuildSlopes(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, double time, TracerStickerNames const& tracerstickersnames) 
{
	const size_t CellNumber = tess.GetPointNo();
	// Get ghost points
	boost::container::flat_map<size_t, ComputationalCell3D> ghost_cells;
	ghost_.operator()(tess, cells, time, tracerstickersnames, ghost_cells);
	// Copy ghost data into new cells vector
	vector<ComputationalCell3D> new_cells(cells);
	new_cells.resize(tess.GetTotalPointNumber());
	for (boost::container::flat_map<size_t, ComputationalCell3D>::const_iterator it = ghost_cells.begin(); it !=
		ghost_cells.end(); ++it)
		new_cells[it->first] = it->second;
	if (SR_)
	{
		size_t Nall = new_cells.size();
		for (size_t j = 0; j < Nall; ++j)
		{
			double gamma = 1.0 / std::sqrt(1 - ScalarProd(new_cells[j].velocity, new_cells[j].velocity));
			new_cells[j].velocity *= gamma;
		}
	}
	// Prepare slopes
	rslopes_.resize(CellNumber, Slope3D(cells[0], cells[0], cells[0]));
	naive_rslopes_.resize(CellNumber);
	Slope3D temp1(cells[0], cells[0], cells[0]);
	ComputationalCell3D temp2(cells[0]);
	ComputationalCell3D temp3(cells[0]);
	ComputationalCell3D temp4(cells[0]);
	ComputationalCell3D temp5(cells[0]);
	vector<ComputationalCell3D> neighbor_list;
	vector<Vector3D> neighbor_mesh_list;
	vector<Vector3D> neighbor_cm_list;
	std::vector<Vector3D> c_ij;
	for (size_t i = 0; i < CellNumber; ++i)
	{
		calc_slope(tess, new_cells, i, slf_, shockratio_, diffusecoeff_, pressure_ratio_, eos_,
			calc_tracers_, naive_rslopes_[i], rslopes_[i], temp1, temp2, temp3, temp4, temp5,
			neighbor_mesh_list, neighbor_cm_list, tracerstickersnames, skip_key_, c_ij, neighbor_list);
	}
#ifdef RICH_MPI
	// communicate ghost slopes
	exchange_ghost_slopes(tess, rslopes_);
#endif //RICH_MPI
}

void LinearGauss3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, double time,
	vector<pair<ComputationalCell3D, ComputationalCell3D> > &res, TracerStickerNames const& tracerstickersnames) const
{
	const size_t CellNumber = tess.GetPointNo();
	vector<size_t> boundaryedges;
	boundaryedges.reserve(static_cast<size_t>(std::pow(static_cast<double>(CellNumber), 0.6666)*8.0));
	// Get ghost points
	boost::container::flat_map<size_t, ComputationalCell3D> ghost_cells;
	ghost_.operator()(tess, cells, time, tracerstickersnames, ghost_cells);
	// Copy ghost data into new cells vector
	vector<ComputationalCell3D> new_cells(cells);
	new_cells.resize(tess.GetTotalPointNumber());
	for (boost::container::flat_map<size_t, ComputationalCell3D>::const_iterator it = ghost_cells.begin(); it !=
		ghost_cells.end(); ++it)
		new_cells[it->first] = it->second;
	if (SR_)
	{
		size_t Nall = new_cells.size();
		for (size_t j = 0; j < Nall; ++j)
		{
			double gamma = 1.0 / std::sqrt(1 - ScalarProd(new_cells[j].velocity, new_cells[j].velocity));
			new_cells[j].velocity *= gamma;
		}
	}
	// Prepare slopes
	rslopes_.resize(CellNumber, Slope3D(cells[0], cells[0], cells[0]));
	naive_rslopes_.resize(CellNumber);
	Slope3D temp1(cells[0], cells[0], cells[0]);
	ComputationalCell3D temp2(cells[0]);
	ComputationalCell3D temp3(cells[0]);
	ComputationalCell3D temp4(cells[0]);
	ComputationalCell3D temp5(cells[0]);
	vector<ComputationalCell3D> neighbor_list;
	vector<Vector3D> neighbor_mesh_list;
	vector<Vector3D> neighbor_cm_list;
	std::vector<Vector3D> c_ij;
	res.resize(tess.GetTotalFacesNumber(), pair<ComputationalCell3D, ComputationalCell3D>(cells[0], cells[0]));
	ComputationalCell3D* cell_ref = 0;
	size_t energy_index = tracerstickersnames.tracer_names.size();
	vector<string>::const_iterator it = binary_find(tracerstickersnames.tracer_names.begin(),
		tracerstickersnames.tracer_names.end(), string("Energy"));
	if (it != tracerstickersnames.tracer_names.end())
		energy_index = static_cast<size_t>(it - tracerstickersnames.tracer_names.begin());
	bool energy_fix = energy_index < tracerstickersnames.tracer_names.size();
	for (size_t i = 0; i < CellNumber; ++i)
	{
		calc_slope(tess, new_cells, i, slf_, shockratio_, diffusecoeff_, pressure_ratio_, eos_,
			calc_tracers_, naive_rslopes_[i], rslopes_[i], temp1, temp2, temp3, temp4, temp5,
			neighbor_mesh_list, neighbor_cm_list, tracerstickersnames, skip_key_, c_ij, neighbor_list);
		face_vec const& faces = tess.GetCellFaces(i);
		const size_t nloop = faces.size();
		for (size_t j = 0; j < nloop; ++j)
		{
			if (tess.GetFaceNeighbors(faces[j]).first == i)
			{
				cell_ref = &res[faces[j]].first;
				ReplaceComputationalCell(*cell_ref, new_cells[i]);
				try
				{
					if(pressure_calc_)
						interp23D(*cell_ref, rslopes_[i], tess.FaceCM(faces[j]), tess.GetCellCM(i), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, rslopes_[i], tess.FaceCM(faces[j]), tess.GetCellCM(i), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
					CheckCell(*cell_ref);
				}
				catch (UniversalError &eo)
				{
					eo.AddEntry("Old density", new_cells[i].density);
					eo.AddEntry("Old internal energy", new_cells[i].internal_energy);
					eo.AddEntry("Face", static_cast<double>(faces[j]));
					eo.AddEntry("Cell", static_cast<double>(i));
					eo.AddEntry("Vx", new_cells[i].velocity.x);
					eo.AddEntry("Vy", new_cells[i].velocity.y);
					eo.AddEntry("Vz", new_cells[i].velocity.z);
					eo.AddEntry("Cell id", static_cast<double>(new_cells[i].ID));
					eo.AddEntry("Interpolated density",cell_ref->density);
					eo.AddEntry("Interpolated pressure",cell_ref->pressure);
					eo.AddEntry("Interpolated internal energy",cell_ref->internal_energy);
					eo.AddEntry("Interpolated Vx",cell_ref->velocity.x);
					eo.AddEntry("Interpolated Vy",cell_ref->velocity.y);
					eo.AddEntry("Interpolated Vz",cell_ref->velocity.z);
					throw eo;
				}
				if (tess.GetFaceNeighbors(faces[j]).second > CellNumber)
					boundaryedges.push_back(faces[j]);
			}
			else
			{
				cell_ref = &res[faces[j]].second;
				ReplaceComputationalCell(*cell_ref, new_cells[i]);
				try
				{
					if (pressure_calc_)
						interp23D(*cell_ref, rslopes_[i], tess.FaceCM(faces[j]), tess.GetCellCM(i), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, rslopes_[i], tess.FaceCM(faces[j]), tess.GetCellCM(i), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
					CheckCell(*cell_ref);
				}
				catch (UniversalError &eo)
				{
					eo.AddEntry("Old density", new_cells[i].density);
					eo.AddEntry("Old internal energy", new_cells[i].internal_energy);
					eo.AddEntry("Face", static_cast<double>(faces[j]));
					eo.AddEntry("Cell", static_cast<double>(i));
					eo.AddEntry("Vx", new_cells[i].velocity.x);
					eo.AddEntry("Vy", new_cells[i].velocity.y);
					eo.AddEntry("Vz", new_cells[i].velocity.z);
					eo.AddEntry("Cell id", static_cast<double>(new_cells[i].ID));
					eo.AddEntry("Interpolated density",cell_ref->density);
					eo.AddEntry("Interpolated pressure",cell_ref->pressure);
					eo.AddEntry("Interpolated internal energy",cell_ref->internal_energy);
					eo.AddEntry("Interpolated Vx",cell_ref->velocity.x);
					eo.AddEntry("Interpolated Vy",cell_ref->velocity.y);
					eo.AddEntry("Interpolated Vz",cell_ref->velocity.z);
					throw eo;
				}
				if (tess.GetFaceNeighbors(faces[j]).first > CellNumber)
					boundaryedges.push_back(faces[j]);
			}
		}
	}
#ifdef RICH_MPI
	// communicate ghost slopes
	exchange_ghost_slopes(tess, rslopes_);
#endif //RICH_MPI
	// Interpolate the boundary edges
	size_t Nboundary = boundaryedges.size();
	for (size_t i = 0; i < Nboundary; ++i)
	{
		size_t N0 = tess.GetFaceNeighbors(boundaryedges[i]).first;
		if (N0 > CellNumber)
		{
			cell_ref = &res[boundaryedges[i]].first;
			ReplaceComputationalCell(*cell_ref, new_cells[N0]);
			try
			{
#ifdef RICH_MPI
				if (tess.BoundaryFace(boundaryedges[i]))
				{ 
					if (pressure_calc_)
						interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
							tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
							tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
				}
				else
				{
					if (pressure_calc_)
						interp23D(*cell_ref, rslopes_[N0], tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, rslopes_[N0], tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
				}
#else
				if (pressure_calc_)
					interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
						tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
				else
				{
					interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
						tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
					if (energy_fix)
						cell_ref->tracers[energy_index] = cell_ref->internal_energy;
				}
#endif //RICH_MPI

				CheckCell(*cell_ref);
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("old density", new_cells[N0].density);
				eo.AddEntry("old internal energy", new_cells[N0].internal_energy);
				eo.AddEntry("Boundary Face", static_cast<double>(boundaryedges[i]));
				eo.AddEntry("Cell", static_cast<double>(N0));
				eo.AddEntry("Vx", new_cells[N0].velocity.x);
				eo.AddEntry("Vy", new_cells[N0].velocity.y);
				eo.AddEntry("Vz", new_cells[N0].velocity.z);
				eo.AddEntry("Cell id", static_cast<double>(new_cells[N0].ID));
				eo.AddEntry("Interpolated density",cell_ref->density);
				eo.AddEntry("Interpolated pressure",cell_ref->pressure);
				eo.AddEntry("Interpolated internal energy",cell_ref->internal_energy);
				eo.AddEntry("Interpolated Vx",cell_ref->velocity.x);
				eo.AddEntry("Interpolated Vy",cell_ref->velocity.y);
				eo.AddEntry("Interpolated Vz",cell_ref->velocity.z);
				eo.AddEntry("Face CMx", tess.FaceCM(boundaryedges[i]).x);
				eo.AddEntry("Face CMy", tess.FaceCM(boundaryedges[i]).y);
				eo.AddEntry("Face CMz", tess.FaceCM(boundaryedges[i]).z);
				eo.AddEntry("Cell CMx", tess.GetCellCM(N0).x);
				eo.AddEntry("Cell CMy", tess.GetCellCM(N0).y);
				eo.AddEntry("Cell CMz", tess.GetCellCM(N0).z);
				size_t N1 = tess.GetFaceNeighbors(boundaryedges[i]).second;
				eo.AddEntry("Other cell ID", static_cast<double>(new_cells[N1].ID));
				eo.AddEntry("Other Cell CMx", tess.GetCellCM(N1).x);
				eo.AddEntry("Other Cell CMy", tess.GetCellCM(N1).y);
				eo.AddEntry("Other Cell CMz", tess.GetCellCM(N1).z);
				eo.AddEntry("Other Cell density", new_cells[N1].density);
				eo.AddEntry("Other Cell pressure", new_cells[N1].pressure);
				eo.AddEntry("Slopex", rslopes_[N0].xderivative.density);
				eo.AddEntry("Slopey", rslopes_[N0].yderivative.density);
				eo.AddEntry("Slopez", rslopes_[N0].zderivative.density);
#ifdef RICH_MPI
				int rank = 0;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				eo.AddEntry("Rank", static_cast<double>(rank));
				for (size_t j = 0; j < tess.GetGhostIndeces().size(); ++j)
					for (size_t k = 0; k < tess.GetGhostIndeces()[j].size(); ++k)
						if (tess.GetGhostIndeces()[j][k] == N0)
						{
							eo.AddEntry("Point recv from proc", static_cast<double>(tess.GetDuplicatedProcs()[j]));
							eo.AddEntry("Point recv index", static_cast<double>(k));
							eo.AddEntry("Point proc index", static_cast<double>(j));
						}
#endif
				throw eo;
			}
		}
		else
		{
			N0 = tess.GetFaceNeighbors(boundaryedges[i]).second;
			cell_ref = &res[boundaryedges[i]].second;
			ReplaceComputationalCell(*cell_ref, new_cells[N0]);
			try
			{
#ifdef RICH_MPI
				if (tess.BoundaryFace(boundaryedges[i]))
				{
					if (pressure_calc_)
						interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
							tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
							tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
				}
				else
				{
					if (pressure_calc_)
						interp23D(*cell_ref, rslopes_[N0], tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
					else
					{
						interp23D(*cell_ref, rslopes_[N0], tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
						if (energy_fix)
							cell_ref->tracers[energy_index] = cell_ref->internal_energy;
					}
				}
#else
				if (pressure_calc_)
					interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
						tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, true);
				else
				{
					interp23D(*cell_ref, ghost_.GetGhostGradient(tess, cells, rslopes_, N0, time, boundaryedges[i],
						tracerstickersnames), tess.FaceCM(boundaryedges[i]), tess.GetCellCM(N0), eos_, tsn_, false);
					if (energy_fix)
						cell_ref->tracers[energy_index] = cell_ref->internal_energy;
				}
#endif //RICH_MPI

				CheckCell(*cell_ref);
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("old density", new_cells[N0].density);
				eo.AddEntry("old internal energy", new_cells[N0].internal_energy);
				eo.AddEntry("Boundary Face", static_cast<double>(boundaryedges[i]));
				eo.AddEntry("Cell", static_cast<double>(N0));
				eo.AddEntry("Vx", new_cells[N0].velocity.x);
				eo.AddEntry("Vy", new_cells[N0].velocity.y);
				eo.AddEntry("Vz", new_cells[N0].velocity.z);
				eo.AddEntry("Cell id", static_cast<double>(new_cells[N0].ID));
				eo.AddEntry("Interpolated density",cell_ref->density);
				eo.AddEntry("Interpolated pressure",cell_ref->pressure);
				eo.AddEntry("Interpolated internal energy",cell_ref->internal_energy);
				eo.AddEntry("Interpolated Vx",cell_ref->velocity.x);
				eo.AddEntry("Interpolated Vy",cell_ref->velocity.y);
				eo.AddEntry("Interpolated Vz",cell_ref->velocity.z);
				eo.AddEntry("Face CMx", tess.FaceCM(boundaryedges[i]).x);
				eo.AddEntry("Face CMy", tess.FaceCM(boundaryedges[i]).y);
				eo.AddEntry("Face CMz", tess.FaceCM(boundaryedges[i]).z);
				eo.AddEntry("Cell CMx", tess.GetCellCM(N0).x);
				eo.AddEntry("Cell CMy", tess.GetCellCM(N0).y);
				eo.AddEntry("Cell CMz", tess.GetCellCM(N0).z);
				size_t N1 = tess.GetFaceNeighbors(boundaryedges[i]).first;
				eo.AddEntry("Other cell ID", static_cast<double>(new_cells[N1].ID));
				eo.AddEntry("Other Cell CMx", tess.GetCellCM(N1).x);
				eo.AddEntry("Other Cell CMy", tess.GetCellCM(N1).y);
				eo.AddEntry("Other Cell CMz", tess.GetCellCM(N1).z);
				eo.AddEntry("Other Cell density", new_cells[N1].density);
				eo.AddEntry("Other Cell pressure", new_cells[N1].pressure);
#ifdef RICH_MPI
				int rank = 0;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				eo.AddEntry("Rank", static_cast<double>(rank));
				for (size_t j = 0; j < tess.GetGhostIndeces().size(); ++j)
					for (size_t k = 0; k < tess.GetGhostIndeces()[j].size(); ++k)
						if (tess.GetGhostIndeces()[j][k] == N0)
						{
							eo.AddEntry("Point recv from proc", static_cast<double>(tess.GetDuplicatedProcs()[j]));
							eo.AddEntry("Point recv index", static_cast<double>(k));
							eo.AddEntry("Point proc index", static_cast<double>(j));
						}
#endif
				throw eo;
			}
		}
	}
	//In SR convert back to velocities
	if (SR_)
	{
		size_t N = res.size();
		for (size_t i = 0; i < N; ++i)
		{
			double factor = 1.0 / std::sqrt(1 + ScalarProd(res[i].first.velocity, res[i].first.velocity));
			res[i].first.velocity *= factor;
			factor = 1.0 / std::sqrt(1 + ScalarProd(res[i].second.velocity, res[i].second.velocity));
			res[i].second.velocity *= factor;
		}
	}
}


vector<Slope3D>& LinearGauss3D::GetSlopes(void)
{
	return rslopes_;
}

vector<Slope3D>& LinearGauss3D::GetSlopesUnlimited(void)const
{
	return naive_rslopes_;
}

