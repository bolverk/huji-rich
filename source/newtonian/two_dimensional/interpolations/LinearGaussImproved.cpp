#include "LinearGaussImproved.hpp"
#include "../../../misc/utils.hpp"
#ifdef RICH_MPI
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#endif //RICH_MPI

typedef pair<ComputationalCell, ComputationalCell> Slope;

namespace 
{
	void GetNeighborMesh(Tessellation const& tess, vector<Edge> const& edges,size_t cell_index,
		vector<Vector2D> &res)
	{
		res.resize(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
				res[i] = tess.GetMeshPoint(neigh1);
			else
				res[i] = tess.GetMeshPoint(neigh0);
		}
	}

	void GetNeighborCM(Tessellation const& tess, vector<Edge> const& edges, size_t cell_index,
		vector<Vector2D> &res)
	{
		res.resize(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			const int neigh0 = edges[i].neighbors.first;
			const int neigh1 = edges[i].neighbors.second;
			if (neigh0 == static_cast<int>(cell_index))
					res[i] = tess.GetCellCM(neigh1);
			else
					res[i] = tess.GetCellCM(neigh0);
		}
	}

	vector<ComputationalCell const*> GetNeighborCells(vector<Edge> const& edges, size_t cell_index,
						   vector<ComputationalCell> const& cells, boost::container::flat_map<size_t, ComputationalCell> const& ghost_cells,size_t /*npoints*/)
	{
		vector<ComputationalCell const*> res(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			size_t other_cell = (edges[i].neighbors.first == static_cast<int>(cell_index)) ? static_cast<size_t>
				(edges[i].neighbors.second) : static_cast<size_t> (edges[i].neighbors.first);
			const boost::container::flat_map<size_t,ComputationalCell>::const_iterator it = 
			  ghost_cells.find(other_cell);
			if(it==ghost_cells.end())
			  res[i]=&cells.at(other_cell);
			else
				res[i]=&it->second;
		}
		return res;
	}

	void GetEdgeList(Tessellation const& tess,
		vector<int> const& edge_indices,vector<Edge> &res)
	{
		res.resize(edge_indices.size());
		for (size_t i = 0; i<edge_indices.size(); ++i)
		{
			res[i] = tess.GetEdge(edge_indices[i]);
		}
	}

	void calc_naive_slope(ComputationalCell const& cell,
		Vector2D const& center, Vector2D const& cell_cm, double cell_volume, vector<ComputationalCell const*> const& neighbors,
		vector<Vector2D> const& neighbor_centers, vector<Vector2D> const& neigh_cm, vector<Edge> const& edge_list,
		pair<ComputationalCell, ComputationalCell> &res)
	{
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
		pair<ComputationalCell, ComputationalCell> vec_compare;
		for (size_t i = 0; i < edge_list.size(); ++i)
		{
			const Vector2D c_ij = CalcCentroid(edge_list[i]) -0.5*(neigh_cm[i] + cell_cm);
			const double e_length = edge_list[i].GetLength();
			const Vector2D r_ij = normalize(neighbor_centers[i] - center)*e_length;
			m[0] -= c_ij.x*r_ij.x;
			m[1] -= c_ij.y*r_ij.x;
			m[2] -= c_ij.x*r_ij.y;
			m[3] -= c_ij.y*r_ij.y;
			if (i == 0)
			{
				vec_compare.first = cell;
				vec_compare.second = cell;
				vec_compare.first *= r_ij.x*0.5;
				vec_compare.second *= r_ij.y*0.5;
			}
			else
			{
				ComputationalCellAddMult(vec_compare.first, cell, r_ij.x*0.5);
				ComputationalCellAddMult(vec_compare.second, cell, r_ij.y*0.5);
			}
			ComputationalCellAddMult(vec_compare.second, *neighbors[i], r_ij.y*0.5);
			ComputationalCellAddMult(vec_compare.first, *neighbors[i], r_ij.x*0.5);
			/*
			vec_compare.first += cell*r_ij.x*0.5;
			vec_compare.first += neighbors[i]*r_ij.x*0.5;
			vec_compare.second += cell*r_ij.y*0.5;
			vec_compare.second += neighbors[i]*r_ij.y*0.5;
			*/
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
		res = vec_compare;
		res.first *= m_inv[0];
		res.second *= m_inv[3];
		ComputationalCellAddMult(res.first, vec_compare.second, m_inv[1]);
		ComputationalCellAddMult(res.second, vec_compare.first, m_inv[2]);
	}

	double PressureRatio(ComputationalCell cell, vector<ComputationalCell const*> const& neigh)
	{
		double res = 1;
		double p = cell.pressure;
		for (size_t i = 0; i<neigh.size(); ++i)
		{
			if (p>neigh[i]->pressure)
				res = std::min(res, neigh[i]->pressure / p);
			else
				res = std::min(res, p / neigh[i]->pressure);
		}
		return res;
	}

	bool is_shock(pair<ComputationalCell,ComputationalCell> const& naive_slope, double cell_width, double shock_ratio,
		ComputationalCell const& cell, vector<ComputationalCell const*> const& neighbor_list, double pressure_ratio,double cs)
	{
		const bool cond1 =(naive_slope.first.velocity.x + naive_slope.second.velocity.y)*
			cell_width<(-shock_ratio)*cs;
		const bool cond2 = PressureRatio(cell, neighbor_list)<pressure_ratio;
		return cond1 || cond2;
	}
	
	ComputationalCell interp(ComputationalCell const& cell,pair<ComputationalCell, ComputationalCell> const& slope,
		Vector2D const& target,Vector2D const& cm)
	{
		ComputationalCell res(cell);
		ComputationalCellAddMult(res, slope.first, target.x - cm.x);
		ComputationalCellAddMult(res, slope.second, target.y - cm.y);
		return res;
	}

	void interp2(ComputationalCell &res, pair<ComputationalCell, ComputationalCell> const& slope,
		Vector2D const& target, Vector2D const& cm)
	{
		ComputationalCellAddMult(res, slope.first, target.x - cm.x);
		ComputationalCellAddMult(res, slope.second, target.y - cm.y);
	}

	void slope_limit(ComputationalCell const& cell,Vector2D const& cm,
		vector<ComputationalCell const*> const& neighbors, vector<Edge> const& edge_list,
		pair<ComputationalCell, ComputationalCell> &slope)
	{
		ComputationalCell cmax(cell), cmin(cell);
		// Find maximum.minimum neighbor values
		for (size_t i = 0; i < neighbors.size(); ++i)
		{
			ComputationalCell const& cell_temp = *neighbors[i];
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
			  (cmax.tracers.begin()+static_cast<int>(j))->second = std::max((cmax.tracers.begin() + static_cast<int>(j))->second,
											(cell_temp.tracers.begin() + static_cast<int>(j))->second);
			  (cmin.tracers.begin() + static_cast<int>(j))->second = std::min((cmin.tracers.begin() + static_cast<int>(j))->second,
											  (cell_temp.tracers.begin() + static_cast<int>(j))->second);
			}
		}
		ComputationalCell maxdiff = cmax - cell,mindiff = cmin - cell;
		// limit the slope
		ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(edge_list[0]), cm);
		ComputationalCell dphi = centroid_val - cell;
		vector<double> psi(4 + cell.tracers.size(), 1);
		for (size_t i = 0; i<edge_list.size(); ++i)
		{
			if (i > 0)
			{
				ReplaceComputationalCell(centroid_val, cell);
				interp2(centroid_val, slope, CalcCentroid(edge_list[i]), cm);
				ReplaceComputationalCell(dphi, centroid_val);
				dphi -= cell;
			}
			//ComputationalCell centroid_val = interp(cell, slope, CalcCentroid(edge_list[i]), cm);
			//ComputationalCell dphi = centroid_val - cell;
			/*Vector2D cent = CalcCentroid(edge_list[i]);
			ComputationalCell dphi(slope.first);
			dphi*=(cent.x - cm.x);
			ComputationalCellAddMult(dphi, slope.second, cent.y - cm.y);*/
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
			for (size_t j = 0; j < dphi.tracers.size(); ++j)
			{
			  double cell_tracer = (cell.tracers.begin() + static_cast<int>(j))->second;
			  double diff_tracer = (maxdiff.tracers.begin() + static_cast<int>(j))->second;
			  if (std::abs((dphi.tracers.begin() + static_cast<int>(j))->second) > 0.1*std::max(std::abs(diff_tracer), std::abs((mindiff.tracers.begin() + static_cast<int>(j))->second)) || (centroid_val.tracers.begin() + static_cast<int>(j))->second *cell_tracer < 0)
				{
				  if ((dphi.tracers.begin() + static_cast<int>(j))->second > std::abs(1e-9*cell_tracer))
				    psi[4 + j] = std::min(psi[4 + j], diff_tracer / (dphi.tracers.begin() + static_cast<int>(j))->second);
					else
					  if ((dphi.tracers.begin() + static_cast<int>(j))->second < -std::abs(1e-9 * cell_tracer))
					    psi[4 + j] = std::min(psi[4 + j], (mindiff.tracers.begin() + static_cast<int>(j))->second
								  / (dphi.tracers.begin() + static_cast<int>(j))->second);
				}
			}
		}
		slope.first.density *= psi[0];
		slope.second.density *= psi[0];
		slope.first.pressure *= psi[1];
		slope.second.pressure *= psi[1];
		slope.first.velocity.x *= psi[2];
		slope.second.velocity.x *= psi[2];
		slope.first.velocity.y *= psi[3];
		slope.second.velocity.y *= psi[3];
		size_t counter = 0;
		for (boost::container::flat_map<std::string, double>::iterator it = slope.first.tracers.begin(); it != slope.first.tracers.end(); ++it)
		{
			it->second *= psi[4 + counter];
	//		safe_retrieve(slope.second.tracers,it->first) *= psi[4 + counter];
			(slope.second.tracers.begin()+static_cast<int>(counter))->second*= psi[4 + counter];
			++counter;
		}
	}

	void shocked_slope_limit(ComputationalCell const& cell, Vector2D const& cm,
		vector<ComputationalCell const*> const& neighbors, vector<Edge> const& edge_list,
		pair<ComputationalCell, ComputationalCell>  &slope,double diffusecoeff)
	{
		ComputationalCell cmax(cell), cmin(cell);
		// Find maximum values
		for (size_t i = 0; i < neighbors.size(); ++i)
		{
			ComputationalCell const& cell_temp = *neighbors[i];
			cmax.density = std::max(cmax.density, cell_temp.density);
			cmax.pressure = std::max(cmax.pressure, cell_temp.pressure);
			cmax.velocity.x = std::max(cmax.velocity.x, cell_temp.velocity.x);
			cmax.velocity.y = std::max(cmax.velocity.y, cell_temp.velocity.y);
			cmin.density = std::min(cmin.density, cell_temp.density);
			cmin.pressure = std::min(cmin.pressure, cell_temp.pressure);
			cmin.velocity.x = std::min(cmin.velocity.x, cell_temp.velocity.x);
			cmin.velocity.y = std::min(cmin.velocity.y, cell_temp.velocity.y);
			for (boost::container::flat_map<std::string, double>::iterator it = cmax.tracers.begin(); it != cmax.tracers.end(); ++it)
			  it->second = std::max(it->second, safe_retrieve(cell_temp.tracers,it->first));
			for (boost::container::flat_map<std::string, double>::iterator it = cmin.tracers.begin(); it != cmin.tracers.end(); ++it)
			  it->second = std::min(it->second, safe_retrieve(cell_temp.tracers,it->first));
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
					psi[0] = std::min(psi[0], std::max(diffusecoeff*(neighbors[i]->density-cell.density)/dphi.density,0.0));
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
			for (boost::container::flat_map<std::string, double>::iterator it = dphi.tracers.begin(); it != dphi.tracers.end(); ++it)
			{
			  double cell_tracer = safe_retrieve(cell.tracers,it->first);
			  double diff_tracer = safe_retrieve(maxdiff.tracers,it->first);
			  double centroid_tracer = safe_retrieve(centroid_val.tracers,it->first);
				if (std::abs(it->second) > 0.1*std::max(std::abs(diff_tracer), std::abs(safe_retrieve(mindiff.tracers,it->first))) || centroid_tracer*cell_tracer < 0)
				{
					if (std::abs(it->second) > std::abs(1e-9*cell_tracer))
					  psi[4 + counter] = std::min(psi[4 + counter], std::max(diffusecoeff*(safe_retrieve(neighbors[i]->tracers,it->first)- cell_tracer) / it->second, 0.0));
				}
				++counter;
			}
		}
		slope.first.density *= psi[0];
		slope.second.density *= psi[0];
		slope.first.pressure *= psi[1];
		slope.second.pressure *= psi[1];
		slope.first.velocity.x *= psi[2];
		slope.second.velocity.x *= psi[2];
		slope.first.velocity.y *= psi[3];
		slope.second.velocity.y *= psi[3];
		size_t counter = 0;
		for (boost::container::flat_map<std::string, double>::iterator it = slope.first.tracers.begin(); it != slope.first.tracers.end(); ++it)
		{
			it->second *= psi[4 + counter];
			safe_retrieve(slope.second.tracers,it->first) *= psi[4 + counter];
			++counter;
		}
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
	 boost::container::flat_map<size_t, ComputationalCell> const& ghost_cells,
	 const vector<string>& flat_tracers,
	 std::pair<ComputationalCell,ComputationalCell> &naive_slope_,
	 std::pair<ComputationalCell, ComputationalCell> & res)
{
	vector<int> edge_indices = tess.GetCellEdges(static_cast<int>(cell_index));
	vector<Edge> edge_list;
	GetEdgeList(tess, edge_indices,edge_list);
	vector<Vector2D> neighbor_mesh_list;
	GetNeighborMesh(tess, edge_list, cell_index,neighbor_mesh_list);
	vector<Vector2D> neighbor_cm_list;
	GetNeighborCM(tess, edge_list, cell_index,neighbor_cm_list);
	vector<ComputationalCell const* > neighbor_list = GetNeighborCells(edge_list, cell_index,cells,ghost_cells,static_cast<size_t>(
		tess.GetPointNo()));

	ComputationalCell const& cell = cells[cell_index];
	calc_naive_slope(cell, tess.GetMeshPoint(static_cast<int>(cell_index)), tess.GetCellCM(static_cast<int>(cell_index)),
		tess.GetVolume(static_cast<int>(cell_index)), neighbor_list, neighbor_mesh_list, neighbor_cm_list, edge_list,res);

	naive_slope_ = res;

	for(size_t i=0;i<flat_tracers.size();++i)
	{
	  res.first.tracers[flat_tracers[i]] = 0;
	  res.second.tracers[flat_tracers[i]] = 0;
	}

	if (slf)
	{
		if (!is_shock(res, tess.GetWidth(static_cast<int>(cell_index)), shockratio, cell, neighbor_list, pressure_ratio,
			      eos.dp2c(cell.density,cell.pressure,cell.tracers)))
		{
			slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, res);
		}
		else
		{
			shocked_slope_limit(cell, tess.GetCellCM(static_cast<int>(cell_index)), neighbor_list, edge_list, res, diffusecoeff);
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
 const vector<string>& flat_tracers): 
  eos_(eos), 
  ghost_(ghost),
  rslopes_(),
  naive_rslopes_(),
  slf_(slf),
  shockratio_(delta_v),
  diffusecoeff_(theta),
  pressure_ratio_(delta_P),
  flat_tracers_(flat_tracers) {}

#ifdef RICH_MPI
namespace
{
	void exchange_ghost_slopes(Tessellation const& tess,vector<Slope > & slopes)
	{
		const boost::mpi::communicator world;
		vector<boost::mpi::request> requests;
		vector<int> const& proc = tess.GetDuplicatedProcs();
		vector<vector<int> > const& tosend = tess.GetDuplicatedPoints();
		vector<vector<Slope> > incoming(proc.size());
		for (size_t i = 0; i < proc.size(); ++i)
		{
			const int dest = proc.at(i);
			requests.push_back(world.isend(dest, 0, VectorValues(slopes, tosend[i])));
			requests.push_back(world.irecv(dest, 0, incoming[i]));
		}
		boost::mpi::wait_all(requests.begin(), requests.end());
		vector<vector<int> > const& ghost_indeces = tess.GetGhostIndeces();
		slopes.resize(tess.GetTotalPointNumber());
		for (size_t i = 0; i < proc.size(); ++i)
		{
			for (size_t j = 0; j < ghost_indeces.at(i).size(); ++j)
				slopes[static_cast<size_t>(ghost_indeces[i][j])] = incoming.at(i).at(j);
		}
	}
}
#endif//RICH_MPI

vector<pair<ComputationalCell, ComputationalCell> > LinearGaussImproved::operator() (const Tessellation& tess,
	const vector<ComputationalCell>& cells,double time) const
{
	const size_t CellNumber = static_cast<size_t>(tess.GetPointNo());
	// Get ghost points
	boost::container::flat_map<size_t,ComputationalCell> ghost_cells = ghost_.operator()(tess,cells,time);
	// Prepare slopes
	rslopes_.resize(CellNumber);
	naive_rslopes_.resize(CellNumber);
	for (size_t i = 0; i<CellNumber; ++i)
	  calc_slope(tess, cells,i,slf_,shockratio_, diffusecoeff_, pressure_ratio_,eos_,ghost_cells,
		flat_tracers_,naive_rslopes_[i],rslopes_[i]);
#ifdef RICH_MPI
	// communicate ghost slopes
	exchange_ghost_slopes(tess, rslopes_);
#endif //RICH_MPI
	// Interpolate the edges
	vector<pair<ComputationalCell, ComputationalCell> > res;
	const size_t edge_number = static_cast<size_t>(tess.GetTotalSidesNumber());
	res.reserve(edge_number);
	for (size_t i = 0; i < edge_number; ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (edge.neighbors.first >= 0 && edge.neighbors.first < static_cast<int>(CellNumber))
		{
			if (edge.neighbors.second >= 0 && edge.neighbors.second < static_cast<int>(CellNumber))
				res.push_back(pair<ComputationalCell, ComputationalCell>(cells[static_cast<size_t>(edge.neighbors.first)],
					cells[static_cast<size_t>(edge.neighbors.second)]));
			else
#ifdef RICH_MPI
			{
				if(tess.GetOriginalIndex(edge.neighbors.second)!= tess.GetOriginalIndex(edge.neighbors.first))
					res.push_back(pair<ComputationalCell, ComputationalCell>(cells[static_cast<size_t>(edge.neighbors.first)],
						cells[static_cast<size_t>(edge.neighbors.second)]));
				else
					res.push_back(pair<ComputationalCell, ComputationalCell>(cells[static_cast<size_t>(edge.neighbors.first)],
						safe_retrieve(ghost_cells, static_cast<size_t>(edge.neighbors.second))));
			}				
#else
				res.push_back(pair<ComputationalCell, ComputationalCell>(cells[static_cast<size_t>(edge.neighbors.first)],
					safe_retrieve(ghost_cells, static_cast<size_t>(edge.neighbors.second))));
#endif
		}
		else
		{
			if (edge.neighbors.second >= 0 && edge.neighbors.second < static_cast<int>(CellNumber))
			{
#ifdef RICH_MPI
				if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
					res.push_back(pair<ComputationalCell, ComputationalCell>(cells[static_cast<size_t>(edge.neighbors.first)],
						cells[static_cast<size_t>(edge.neighbors.second)]));
				else
					res.push_back(pair<ComputationalCell, ComputationalCell>(safe_retrieve(ghost_cells, static_cast<size_t>(edge.neighbors.first)),
						cells[static_cast<size_t>(edge.neighbors.second)]));
#else
				res.push_back(pair<ComputationalCell, ComputationalCell>(safe_retrieve(ghost_cells, static_cast<size_t>(edge.neighbors.first)),
					cells[static_cast<size_t>(edge.neighbors.second)]));
#endif
			}
			else
				throw UniversalError("Both sides of edge are ghost");
		}		
	}
	for (size_t i = 0; i < edge_number; ++i)
	{
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (edge.neighbors.first >= 0 && edge.neighbors.first < static_cast<int>(CellNumber)
#ifdef RICH_MPI
			|| tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second)
#endif
			)
		{
			interp2(res[i].first, rslopes_[static_cast<size_t>(edge.neighbors.first)],
				CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
		}
		else
		{
		  interp2(res[i].first, ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
			  edge.neighbors.first), time,edge), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.first));
		}
		if (edge.neighbors.second >= 0 && edge.neighbors.second < static_cast<int>(CellNumber)
#ifdef RICH_MPI
			|| tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second)
#endif
			)
		{
			interp2(res[i].second, rslopes_[static_cast<size_t>(edge.neighbors.second)],
				CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
		}
		else
		{
		  const ComputationalCell& cell = safe_retrieve
		    (ghost_cells,
		     static_cast<size_t>(edge.neighbors.second));
			res[i].second = interp(cell, ghost_.GetGhostGradient(tess, cells, rslopes_, static_cast<size_t>(
				edge.neighbors.second),time,edge), CalcCentroid(edge), tess.GetCellCM(edge.neighbors.second));
		}
	}
	return res;
}


vector<pair<ComputationalCell, ComputationalCell> >& LinearGaussImproved::GetSlopes(void)const
{
	return rslopes_;
}

vector<pair<ComputationalCell, ComputationalCell> >& LinearGaussImproved::GetSlopesUnlimited(void)const
{
	return naive_rslopes_;
}

