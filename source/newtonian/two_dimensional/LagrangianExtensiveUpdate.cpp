#include "LagrangianExtensiveUpdate.hpp"
#include "../../misc/utils.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif
namespace
{
	bool bracketed(int low, int arg, int high)
	{
		return arg >= low && high > arg;
	}
}

LagrangianExtensiveUpdate::LagrangianExtensiveUpdate(LagrangianFlux const& lflux, GhostPointGenerator
	const& ghost) :lflux_(lflux), ghost_(ghost) {}

void LagrangianExtensiveUpdate::operator()(const vector<Extensive>& fluxes, const PhysicalGeometry& /*pg*/,
	const Tessellation& tess, const double dt, const CacheData& cd, const vector<ComputationalCell>& cells,
	vector<Extensive>& extensives, double time, TracerStickerNames const& tracerstickernames) const
{
	const vector<Edge>& edge_list = tess.getAllEdges();
	Extensive delta = dt*cd.areas[0] * fluxes[0];
	size_t N = static_cast<size_t>(tess.GetPointNo());
	std::vector<double> dA(N, 0), dWs(N, 0);
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	//assert(indexX < tracerstickernames.tracer_names.size() && indexY < tracerstickernames.tracer_names.size());
	bool area_calc = false;
	if (indexX < tracerstickernames.tracer_names.size())
	{
		area_calc = true;
		for (size_t i = 0; i < N; ++i)
		{
			extensives[i].tracers[indexX] *= 0.85;
			extensives[i].tracers[indexY] *= 0.85;
		}
	}

	boost::container::flat_map<size_t,ComputationalCell> ghosts = ghost_(tess, cells, time, tracerstickernames);
	vector<ComputationalCell> newcells(cells);
	newcells.resize(static_cast<size_t>(tess.GetTotalPointNumber()));
	for (boost::container::flat_map<size_t, ComputationalCell>::iterator it = ghosts.begin(); it != ghosts.end(); ++it)
		newcells[it->first] = it->second;

	for (size_t i = 0; i < edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		ReplaceExtensive(delta, fluxes[i]);
		delta *= dt*cd.areas[i];
		double deltaWs = -lflux_.ws_[i] * cd.areas[i] * dt;
		Vector2D normal = normalize(tess.GetMeshPoint(edge.neighbors.second) - tess.GetMeshPoint(edge.neighbors.first));
		if (lflux_.Lag_calc_[i])
		{
			double p_star = ScalarProd(fluxes[i].momentum, normal);
			double v_star = fluxes[i].energy / p_star;
			double v_new = (v_star - lflux_.ws_[i]);
			if (v_new*v_star > 0)
			{
				if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
					&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.second)].pressure)
					v_new = std::max(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.second)].velocity,
						normal));
				if (v_new < 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
					&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.first)].pressure)
					v_new = std::min(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.first)].velocity,
						normal));
				delta.energy = p_star*v_new*cd.areas[i] * dt;
			}
			else
			{
				if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second))
					delta.energy = cd.areas[i] * dt*v_new*newcells[static_cast<size_t>(edge.neighbors.first)].pressure;
				else
					if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
						delta.energy = cd.areas[i] * dt*v_new*newcells[static_cast<size_t>(edge.neighbors.second)].pressure;
					else
						delta.energy = 0;
			}
		}
		if (bracketed(0, edge.neighbors.first, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.first)] -= delta;
			if (area_calc)
			{
				extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexX] -= normal.x*deltaWs;
				extensives[static_cast<size_t>(edge.neighbors.first)].tracers[indexY] -= normal.y*deltaWs;
			}
		}
		if (bracketed(0, edge.neighbors.second, tess.GetPointNo()))
		{
			extensives[static_cast<size_t>(edge.neighbors.second)] += delta;
			if (area_calc)
			{
				extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexX] -= normal.x*deltaWs;
				extensives[static_cast<size_t>(edge.neighbors.second)].tracers[indexY] -= normal.y*deltaWs;
			}
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		if ((!(extensives[i].mass > 0)) || (!(extensives[i].energy > 0)) || (!std::isfinite(extensives[i].momentum.x)) || (!std::isfinite(extensives[i].momentum.y)))
		{
			int rank = 0;
#ifdef RICH_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
			std::cout << "Bad cell in LagrangianExtensiveUpdate, cell " << i << " rank " << rank << " time " << time << std::endl;
			std::cout << "mass " << extensives[i].mass << " energy " << extensives[i].energy << " Entropy " <<
				extensives[i].tracers[0] << " momentum" << abs(extensives[i].momentum) << " volume " << cd.volumes[i]<< std::endl;
			std::cout << "Old cell, density " << cells[i].density << " pressure " << cells[i].pressure << " v " <<
				abs(cells[i].velocity) << std::endl;
			std::vector<int> temp = tess.GetCellEdges(static_cast<int>(i));
			for (size_t j = 0; j < temp.size(); ++j)
			{
				Edge edge = tess.GetEdge(temp[j]);
				int N0 = tess.GetEdge(temp[j]).neighbors.first;
				int N1 = tess.GetEdge(temp[j]).neighbors.second;
				double Area = tess.GetEdge(temp[j]).GetLength() * dt;
				std::cout << "Edge " << temp[j] << " neigh " << N0 << "," << N1 << " mass=" << fluxes[temp[j]].mass*Area <<
					" energy " << fluxes[temp[j]].energy*Area << " momentum=" << abs(fluxes[temp[j]].momentum)*Area <<
					" L*dt " << Area << " N0 " << tess.GetMeshPoint(N0).x << "," << tess.GetMeshPoint(N0).y << " N1 " 
					<< tess.GetMeshPoint(N1).x << "," << tess.GetMeshPoint(N1).y << std::endl;
				std::cout << "dl " << newcells[N0].density << " pl " << newcells[N0].pressure << " vxl " << newcells[N0].velocity.x << " vyl " << newcells[N0].velocity.y  << std::endl;
				std::cout << "dr " << newcells[N1].density << " pr " << newcells[N1].pressure << " vxr " << newcells[N1].velocity.x << " vyr " << newcells[N1].velocity.y  << std::endl;
				Vector2D normal = normalize(tess.GetMeshPoint(N1) - tess.GetMeshPoint(N0));
				std::cout << "dx " << normal.x << " dy " << normal.y << std::endl;
				if (lflux_.Lag_calc_[temp[j]])
				{
					double p_star = ScalarProd(fluxes[temp[j]].momentum, normal);
					double v_star = fluxes[temp[j]].energy / p_star;
					double v_new = (v_star - lflux_.ws_[temp[j]]);
					std::cout << "Old pstar " << p_star << " vstar " << v_star << " vnew " << v_new << std::endl;
					if (v_new*v_star > 0)
					{
						if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
							&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.second)].pressure)
							v_new = std::max(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.second)].velocity,
								normal));
						if (v_new < 0 && tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first)
							&& p_star > 1.2*newcells[static_cast<size_t>(edge.neighbors.first)].pressure)
							v_new = std::min(v_new, ScalarProd(newcells[static_cast<size_t>(edge.neighbors.first)].velocity,
								normal));
					}
					else
					{
						if (v_new > 0 && tess.GetOriginalIndex(edge.neighbors.first) != tess.GetOriginalIndex(edge.neighbors.second))
							p_star = newcells[static_cast<size_t>(edge.neighbors.first)].pressure;
						else
							if (tess.GetOriginalIndex(edge.neighbors.second) != tess.GetOriginalIndex(edge.neighbors.first))
								p_star = newcells[static_cast<size_t>(edge.neighbors.second)].pressure;
							else
								p_star = 0;
					}
					std::cout << "Pstar " << p_star << " Ustar " << v_new << std::endl;
				}
			}
			assert(false);
		}
	}
}

