#include "LinearGaussImproved.hpp"
#include "../../../misc/utils.hpp"

namespace 
{
	vector<Vector2D> GetNeighborMesh(Tessellation const& tess, vector<Edge> const& edges,size_t cell_index)
	{
		vector<Vector2D> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
				res[i] = tess.GetMeshPoint(neigh1);
			else
				res[i] = tess.GetMeshPoint(neigh0);
		}
		return res;
	}

	vector<Vector2D> GetNeighborCM(Tessellation const& tess, vector<Edge> const& edges, size_t cell_index)
	{
		vector<Vector2D> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
					res[i] = tess.GetCellCM(neigh1);
			else
					res[i] = tess.GetCellCM(neigh0);
		}
		return res;
	}

	vector<ComputationalCell> GetNeighborCells(vector<Edge> const& edges, size_t cell_index,
		vector<ComputationalCell> const& cells, boost::container::flat_map<size_t, ComputationalCell> const& ghost_cells,size_t npoints)
	{
		vector<ComputationalCell> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			size_t other_cell = (edges[i].neighbors.first == static_cast<int>(cell_index)) ? static_cast<size_t>
				(edges[i].neighbors.second) : static_cast<size_t> (edges[i].neighbors.first);
			if (other_cell < npoints)
				res[i] = cells[other_cell];
			else
			{
				boost::container::flat_map<size_t, ComputationalCell>::const_iterator it = ghost_cells.find(other_cell);
				if (it == ghost_cells.end())
					throw UniversalError("No matching ghost cell");
				res[i] = it->second;
			}
		}
		return res;
	}

	vector<Edge> GetEdgeList(Tessellation const& tess,
		vector<int> const& edge_indices)
	{
		vector<Edge> res(edge_indices.size());
		for (size_t i = 0; i<edge_indices.size(); ++i)
		{
			res[i] = tess.GetEdge(edge_indices[i]);
		}
		return res;
	}

	pair<ComputationalCell, ComputationalCell> calc_naive_slope(ComputationalCell const& cell,
		Vector2D const& center, Vector2D const& cell_cm, double cell_volume, vector<ComputationalCell> const& neighbors,
		vector<Vector2D> const& neighbor_centers, vector<Vector2D> const& neigh_cm, vector<Edge> const& edge_list)
	{
		pair<ComputationalCell, ComputationalCell> res;
		size_t n = edge_list.size();
		if (n>20)
		{
			UniversalError eo("Cell has too many neighbors");
			eo.AddEntry("Cell x cor", center.x);
			eo.AddEntry("Cell y cor", center.y);
			throw eo;
		}
		// Create the matrix to invert and the vector to compare
		vector<double> m(4, 0);
		pair<ComputationalCell,ComputationalCell> vec_compare
		  (0*cell,0*cell);
		for (size_t i = 0; i < edge_list.size(); ++i)
		{
			const Vector2D c_ij = CalcCentroid(edge_list[i]) -0.5*(neigh_cm[i] + cell_cm);
			const double e_length = edge_list[i].GetLength();
			const Vector2D r_ij = normalize(neighbor_centers[i] - center)*e_length;
			m[0] -= c_ij.x*r_ij.x;
			m[1] -= c_ij.y*r_ij.x;
			m[2] -= c_ij.x*r_ij.y;
			m[3] -= c_ij.y*r_ij.y;
			vec_compare.first = vec_compare.first + (cell + neighbors[i])*r_ij.x*0.5;
			vec_compare.second = vec_compare.second + (cell + neighbors[i])*r_ij.y*0.5;
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
		return pair<ComputationalCell, ComputationalCell>(m_inv[0] * vec_compare.first + m_inv[1] * vec_compare.second,
			m_inv[2] * vec_compare.first + m_inv[3] * vec_compare.second);
	}

	double PressureRatio(ComputationalCell cell, vector<ComputationalCell> const& neigh)
	{
		double res = 1;
		double p = cell.pressure;
		for (size_t i = 0; i<neigh.size(); ++i)
		{
			if (p>neigh[i].pressure)
				res = std::min(res, neigh[i].pressure / p);
			else
				res = std::min(res, p / neigh[i].pressure);
		}
		return res;
	}

	bool is_shock(pair<ComputationalCell,ComputationalCell> const& naive_slope, double cell_width, double shock_ratio,
		ComputationalCell const& cell, vector<ComputationalCell> const& neighbor_list, double pressure_ratio,double cs)
	{
		const bool cond1 =(naive_slope.first.velocity.x + naive_slope.second.velocity.y)*
			cell_width<(-shock_ratio)*cs;
		const bool cond2 = PressureRatio(cell, neighbor_list)<pressure_ratio;
		return cond1 || cond2;
	}
	
	ComputationalCell interp(ComputationalCell const& cell,pair<ComputationalCell, ComputationalCell> const& slope,
		Vector2D const& target,Vector2D const& cm)
	{
		return cell + (target.x - cm.x)*slope.first + (target.y - cm.y)*slope.second;
	}

	pair<ComputationalCell, ComputationalCell> slope_limit(ComputationalCell const& cell,Vector2D const& cm,
		vector<ComputationalCell> const& neighbors, vector<Edge> const& edge_list,
		pair<ComputationalCell, ComputationalCell> const& slope)
	{
		pair<ComputationalCell, ComputationalCell> res(slope);
		ComputationalCell cmax(cell), cmin(cell);
		// Find maximum.minimum neighbor values
		for (size_t i = 0; i < neighbors.size(); ++i)
		{
			ComputationalCell const& cell_temp = neighbors[i];
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			for (boost::container::flat_map<std::string, double>::iterator it = cmax.tracers.begin(); it != cmax.tracers.end(); ++it)
				it->second = std::max(it->second, cell_temp.tracers.find(it->first)->second);
			for (boost::container::flat_map<std::string, double>::iterator it = cmin.tracers.begin(); it != cmin.tracers.end(); ++it)
				it->second = std::min(it->second, cell_temp.tracers.find(it->first)->second);
		}
		ComputationalCell maxdiff = cmax - cell,mindiff = cmin - cell;
		// limit the slope
		vector<double> psi(4 + cell.tracers.size(), 1);
		for (size_t i = 0; i<edge_list.size(); ++i)
		{
			ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(edge_list[i]), cm);
			ComputationalCell dphi = centroid_val - cell;
			// density
			if (std::abs(dphi.density) > 0.1*std::max(std::abs(maxdiff.density),std::abs(mindiff.density))||centroid_val.density*cell.density < 0)
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
			if (std::abs(dphi.velocity.y) > 0.1*std::max(std::abs(maxdiff.velocity.y),std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (dphi.velocity.y > std::abs(1e-9*cell.velocity.y))
					psi[3] = std::min(psi[3], maxdiff.velocity.y / dphi.velocity.y);
				else
					if (dphi.velocity.y<-std::abs(1e-9*cell.velocity.y))
						psi[3] = std::min(psi[3], mindiff.velocity.y / dphi.velocity.y);
			}
			// tracers
			size_t counter = 0;
			for (boost::container::flat_map<std::string, double>::iterator it = dphi.tracers.begin(); it != dphi.tracers.end(); ++it)
			{
				double cell_tracer = cell.tracers.find(it->first)->second;
				double diff_tracer = maxdiff.tracers.find(it->first)->second;
				double centroid_tracer = centroid_val.tracers.find(it->first)->second;
				if (std::abs(it->second) > 0.1*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers.find(it->first)->second)) || centroid_tracer*cell_tracer < 0)
				{
					if (it->second > std::abs(1e-9*cell_tracer))
						psi[4 + counter] = std::min(psi[4 + counter], diff_tracer / it->second);
					else
						if (it->second < -std::abs(1e-9 * cell_tracer))
						  psi[4 + counter] = std::min(psi[4 + counter], safe_retrieve(mindiff.tracers,it->first)
							/ it->second);
				}
				++counter;
			}
		}
		res.first.density *= psi[0];
		res.second.density *= psi[0];
		res.first.pressure *= psi[1];
		res.second.pressure *= psi[1];
		res.first.velocity.x *= psi[2];
		res.second.velocity.x *= psi[2];
		res.first.velocity.y *= psi[3];
		res.second.velocity.y *= psi[3];
		size_t counter = 0;
		for (boost::container::flat_map<std::string, double>::iterator it = res.first.tracers.begin(); it != res.first.tracers.end(); ++it)
		{
			it->second *= psi[4 + counter];
			res.second.tracers.find(it->first)->second *= psi[4 + counter];
			++counter;
		}
		return res;
	}

	pair<ComputationalCell, ComputationalCell> shocked_slope_limit(ComputationalCell const& cell, Vector2D const& cm,
		vector<ComputationalCell> const& neighbors, vector<Edge> const& edge_list,
		pair<ComputationalCell, ComputationalCell> const& slope,double diffusecoeff)
	{
		pair<ComputationalCell, ComputationalCell> res(slope);
		ComputationalCell cmax(cell), cmin(cell);
		// Find maximum values
		for (size_t i = 0; i < neighbors.size(); ++i)
		{
			ComputationalCell const& cell_temp = neighbors[i];
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			for (boost::container::flat_map<std::string, double>::iterator it = cmax.tracers.begin(); it != cmax.tracers.end(); ++it)
				it->second = std::max(it->second, cell_temp.tracers.find(it->first)->second);
			for (boost::container::flat_map<std::string, double>::iterator it = cmin.tracers.begin(); it != cmin.tracers.end(); ++it)
				it->second = std::min(it->second, cell_temp.tracers.find(it->first)->second);
		}
		ComputationalCell maxdiff = cmax - cell, mindiff = cmin - cell;
		// limit the slope
		vector<double> psi(4 + cell.tracers.size(), 1);
		for (size_t i = 0; i<edge_list.size(); ++i)
		{
			ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(edge_list[i]), cm);
			ComputationalCell dphi = centroid_val - cell;
			// density
			if (std::abs(dphi.density) > 0.1*std::max(std::abs(maxdiff.density), std::abs(mindiff.density)) || centroid_val.density*cell.density < 0)
			{
				if (std::abs(dphi.density) > 1e-9*cell.density)
					psi[0] = std::min(psi[0], std::max(diffusecoeff*(neighbors[i].density-cell.density)/dphi.density,0.0));
			}
			// pressure
			if (std::abs(dphi.pressure) > 0.1*std::max(std::abs(maxdiff.pressure), std::abs(mindiff.pressure)) || centroid_val.pressure*cell.pressure < 0)
			{
				if (std::abs(dphi.pressure) > 1e-9*cell.pressure)
					psi[1] = std::min(psi[1], std::max(diffusecoeff*(neighbors[i].pressure - cell.pressure) / dphi.pressure, 0.0));
			}
			// xvelocity
			if (std::abs(dphi.velocity.x) > 0.1*std::max(std::abs(maxdiff.velocity.x), std::abs(mindiff.velocity.x)) || centroid_val.velocity.x*cell.velocity.x < 0)
			{
				if (std::abs(dphi.velocity.x) > 1e-9*cell.velocity.x)
					psi[2] = std::min(psi[2], std::max(diffusecoeff*(neighbors[i].velocity.x - cell.velocity.x) / dphi.velocity.x, 0.0));
			}
			// yvelocity
			if (std::abs(dphi.velocity.y) > 0.1*std::max(std::abs(maxdiff.velocity.y), std::abs(mindiff.velocity.y)) || centroid_val.velocity.y*cell.velocity.y < 0)
			{
				if (std::abs(dphi.velocity.y) > 1e-9*cell.velocity.y)
					psi[3] = std::min(psi[3], std::max(diffusecoeff*(neighbors[i].velocity.y - cell.velocity.y) / dphi.velocity.y, 0.0));
			}
			// tracers
			size_t counter = 0;
			for (boost::container::flat_map<std::string, double>::iterator it = dphi.tracers.begin(); it != dphi.tracers.end(); ++it)
			{
				double cell_tracer = cell.tracers.find(it->first)->second;
				double diff_tracer = maxdiff.tracers.find(it->first)->second;
				double centroid_tracer = centroid_val.tracers.find(it->first)->second;
				if (std::abs(it->second) > 0.1*std::max(std::abs(diff_tracer), std::abs(mindiff.tracers.find(it->first)->second)) || centroid_tracer*cell_tracer < 0)
				{
					if (std::abs(it->second) > std::abs(1e-9*cell_tracer))
						psi[4 + counter] = std::min(psi[4 + counter], std::max(diffusecoeff*(neighbors[i].tracers.find(it->first)->second - cell_tracer) / it->second, 0.0));
				}
				++counter;
			}
		}
		res.first.density *= psi[0];
		res.second.density *= psi[0];
		res.first.pressure *= psi[1];
		res.second.pressure *= psi[1];
		res.first.velocity.x *= psi[2];
		res.second.velocity.x *= psi[2];
		res.first.velocity.y *= psi[3];
		res.second.velocity.y *= psi[3];
		size_t counter = 0;
		for (boost::container::flat_map<std::string, double>::iterator it = res.first.tracers.begin(); it != res.first.tracers.end(); ++it)
		{
			it->second *= psi[4 + counter];
			safe_retrieve(res.second.tracers,it->first) *= psi[4 + counter];
			++counter;
		}
		return res;
	}

	pair<ComputationalCell,ComputationalCell> calc_slope
	(Tessellation const& tess,
	 vector<ComputationalCell> const& cells,
	 size_t cell_index,
	 bool slf,
	 double shockratio,
	 double diffusecoeff,
	 double pressure_ratio,
	 EquationOfState const& eos,
	 boost::container::flat_map<size_t, ComputationalCell> const& ghost_cells,
	 const vector<string>& flat_tracers)
{
	vector<int> edge_indices = tess.GetCellEdges(static_cast<int>(cell_index));
	vector<Edge> edge_list = GetEdgeList(tess, edge_indices);
	vector<Vector2D> neighbor_mesh_list = GetNeighborMesh(tess, edge_list,cell_index);
	vector<Vector2D> neighbor_cm_list = GetNeighborCM(tess, edge_list,cell_index);
	vector<ComputationalCell> neighbor_list = GetNeighborCells(edge_list, cell_index,cells,ghost_cells,static_cast<size_t>(
		tess.GetPointNo()));

	pair<ComputationalCell, ComputationalCell> naive_slope, s_compare;
	ComputationalCell const& cell = cells[cell_index];
	naive_slope = calc_naive_slope(cell, tess.GetMeshPoint(static_cast<int>(cell_index)), tess.GetCellCM(static_cast<int>(cell_index)),
		tess.GetVolume(static_cast<int>(cell_index)), neighbor_list, neighbor_mesh_list, neighbor_cm_list, edge_list);

	for(size_t i=0;i<flat_tracers.size();++i){
	  naive_slope.first.tracers[flat_tracers[i]] = 0;
	  naive_slope.second.tracers[flat_tracers[i]] = 0;
	}

	if (slf)
	{
		if (!is_shock(naive_slope, tess.GetWidth(static_cast<int>(cell_index)), shockratio, cell, neighbor_list, pressure_ratio,
			      eos.dp2c(cell.density,cell.pressure,cell.tracers)))
		{
			return slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, naive_slope);
		}
		else
		{
			return shocked_slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, naive_slope, diffusecoeff);
		}
	}
	else
	{
		return naive_slope;
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
 const vector<string>& flat_tracers): 
  eos_(eos), 
  ghost_(ghost),
  rslopes_(),
  slf_(slf),
  shockratio_(delta_v),
  diffusecoeff_(theta),
  pressure_ratio_(delta_P),
  flat_tracers_(flat_tracers) {}

vector<pair<ComputationalCell, ComputationalCell> > LinearGaussImproved::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells) const
{
	const size_t CellNumber = static_cast<size_t>(tess.GetPointNo());
	// Get ghost points
	boost::container::flat_map<size_t,ComputationalCell> ghost_cells = ghost_.operator()(tess,cells);
	// Prepare slopes
	rslopes_.resize(CellNumber);
	for (size_t i = 0; i<CellNumber; ++i)
	  rslopes_[i] = calc_slope(tess, cells,i,slf_,shockratio_, diffusecoeff_, pressure_ratio_,eos_,ghost_cells, flat_tracers_);
	// Interpolate the edges
	vector<pair<ComputationalCell, ComputationalCell> > res;
	const size_t edge_number = static_cast<size_t>(tess.GetTotalSidesNumber());
	res.resize(edge_number);
	for (size_t i = 0; i < edge_number; ++i)
	{
		pair<ComputationalCell, ComputationalCell> cell_temp;
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (edge.neighbors.first >= 0 && edge.neighbors.first < static_cast<int>(CellNumber))
			cell_temp.first = interp(cells[static_cast<size_t>(edge.neighbors.first)], rslopes_[static_cast<size_t>(edge.neighbors.first)],
			CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
		else
		{
		  ComputationalCell const& cell = safe_retrieve(ghost_cells,static_cast<size_t>(edge.neighbors.first));
			cell_temp.first = interp(cell, ghost_.GetGhostGradient(tess,cells,rslopes_,static_cast<size_t>(
				edge.neighbors.first)),CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
		}
		if (edge.neighbors.second >= 0 && edge.neighbors.second < static_cast<int>(CellNumber))
			cell_temp.second = interp(cells[static_cast<size_t>(edge.neighbors.second)], rslopes_[static_cast<size_t>(edge.neighbors.second)],
			CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
		else
		{
		  const ComputationalCell& cell = safe_retrieve
		    (ghost_cells,
		     static_cast<size_t>(edge.neighbors.second));
			cell_temp.second = interp(cell, ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
				edge.neighbors.second)), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
		}
		res[i] = cell_temp;
	}
	return res;
}


vector<pair<ComputationalCell, ComputationalCell> > LinearGaussImproved::GetSlopes(void)const
{
	return rslopes_;
}
