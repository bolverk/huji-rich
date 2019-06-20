#include "AMR3D.hpp"
#include "../../3D/GeometryCommon/Voronoi3D.hpp"
#include "../../misc/utils.hpp"
#include <boost/array.hpp>
#include <iostream>
#include "../../3D/r3d/Intersection3D.hpp"
#include <boost/scoped_ptr.hpp>

//#define debug_amr 1

namespace
{
	void RemoveRefineNeighborRemove(Tessellation3D const& tess, std::vector<size_t> const& remove,
		std::vector<size_t> &refine, std::vector<Vector3D> &refine_direction)
	{
		vector<size_t> neigh;
		size_t Nrefine = refine.size();
		vector<size_t> newrefine;
		newrefine.reserve(Nrefine);
		std::vector<Vector3D> new_direction;
		new_direction.reserve(Nrefine);
		for (size_t i = 0; i < Nrefine; ++i)
		{
			tess.GetNeighbors(refine[i], neigh);
			neigh.push_back(refine[i]);
			bool good = true;
			for (size_t j = 0; j < neigh.size(); ++j)
			{
				if (std::binary_search(remove.begin(), remove.end(), neigh[j]))
				{
					good = false;
					break;
				}
			}
			if (good)
			{
				newrefine.push_back(refine[i]);
				if (!refine_direction.empty())
					new_direction.push_back(refine_direction[i]);
			}
		}
		refine = newrefine;
		refine_direction = new_direction;
	}

	std::pair<vector<size_t>, vector<double> > RemoveNeighbors(vector<double> const& merits, vector<size_t> const& candidates,
		Tessellation3D const& tess)
	{
		vector<size_t> result_names;
		vector<double> result_merits;
		if (merits.size() != candidates.size())
			throw UniversalError("Merits and Candidates don't have same size in RemoveNeighbors");
		// Make sure there are no neighbors
		vector<size_t> bad_neigh;
		vector<size_t> neigh;
		for (size_t i = 0; i < merits.size(); ++i)
		{
			bool good = true;
			tess.GetNeighbors(candidates[i], neigh);
			size_t nneigh = neigh.size();
			if (find(bad_neigh.begin(), bad_neigh.end(), candidates[i]) != bad_neigh.end())
				good = false;
			else
			{
				for (size_t j = 0; j < nneigh; ++j)
				{
					if (binary_search(candidates.begin(), candidates.end(), neigh[j]))
					{
						if (merits[i] < merits[static_cast<size_t>(lower_bound(candidates.begin(), candidates.end(),
							neigh[j]) - candidates.begin())])
						{
							good = false;
							break;
						}
						if (fabs(merits[i] - merits[static_cast<size_t>(lower_bound(candidates.begin(), candidates.end(),
							neigh[j]) - candidates.begin())]) < 1e-9)
						{
							if (find(bad_neigh.begin(), bad_neigh.end(), neigh[j]) == bad_neigh.end())
								bad_neigh.push_back(neigh[j]);
						}
					}
				}
			}
			if (good)
			{
				result_names.push_back(candidates[i]);
				result_merits.push_back(merits[i]);
			}
		}
		return std::pair<vector<size_t>, vector<double> >(result_names, result_merits);
	}

#ifdef RICH_MPI
	std::pair<vector<size_t>, vector<double> > RemoveMPINeighbors(vector<double> const& merits, vector<size_t> const& candidates,
		Tessellation3D const& tess)
	{
		vector<vector<size_t> > duplicated_indeces = tess.GetDuplicatedPoints();
		vector<vector<size_t> > sort_indeces(duplicated_indeces.size());
		for (size_t i = 0; i < duplicated_indeces.size(); ++i)
		{
			sort_index(duplicated_indeces[i], sort_indeces[i]);
			std::sort(duplicated_indeces[i].begin(), duplicated_indeces[i].end());
		}
		size_t nproc = duplicated_indeces.size();
		vector<vector<size_t> > indeces(nproc);
		vector<vector<double> > merit_send(nproc);
		vector<size_t> neigh;
		size_t Norg = tess.GetPointNo();
		size_t Nreomve = merits.size();
		// Send/recv data
		for (size_t i = 0; i < Nreomve; ++i)
		{
			tess.GetNeighbors(candidates[i], neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] >= Norg)
				{
					for (size_t k = 0; k < nproc; ++k)
					{
						vector<size_t>::const_iterator it = binary_find(duplicated_indeces[k].begin(),
							duplicated_indeces[k].end(), candidates[i]);
						if (it != duplicated_indeces[k].end())
						{
							indeces[k].push_back(sort_indeces[k][static_cast<size_t>(it - duplicated_indeces[k].begin())]);
							merit_send[k].push_back(merits[i]);
						}
					}
				}
			}
		}
		indeces = MPI_exchange_data(tess.GetDuplicatedProcs(), indeces);
		merit_send = MPI_exchange_data(tess.GetDuplicatedProcs(), merit_send);
		vector<size_t> all_indeces, temp;
		vector<double> all_merits;
		for (size_t i = 0; i < nproc; ++i)
		{
			if (!indeces[i].empty())
			{
				indeces[i] = VectorValues(tess.GetGhostIndeces()[i], indeces[i]);
				all_indeces.insert(all_indeces.end(), indeces[i].begin(), indeces[i].end());
				all_merits.insert(all_merits.end(), merit_send[i].begin(), merit_send[i].end());
			}
		}
		sort_index(all_indeces, temp);
		std::sort(all_indeces.begin(), all_indeces.end());
		all_merits = VectorValues(all_merits, temp);
		// remove neighbors
		std::pair<vector<size_t>, vector<double> > res;
		res.first.reserve(Nreomve);
		res.second.reserve(Nreomve);
		for (size_t i = 0; i < Nreomve; ++i)
		{
			bool good = true;
			tess.GetNeighbors(candidates[i], neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] >= Norg)
				{
					vector<size_t>::const_iterator it = binary_find(all_indeces.begin(), all_indeces.end(), neigh[j]);
					if (it != all_indeces.end())
					{
						if (all_merits[static_cast<size_t>(it - all_indeces.begin())] > merits[i])
						{
							good = false;
							break;
						}
					}
				}
			}
			if (good)
			{
				res.first.push_back(candidates[i]);
				res.second.push_back(merits[i]);
			}
		}
		return res;
	}
#endif //RICH_MPI

	std::vector<Vector3D> GetNewPoints(Tessellation3D const& tess, std::pair<vector<size_t>,
		vector<Vector3D> > &ToRefine
#ifdef RICH_MPI
		, Tessellation3D const& tproc
#endif
	)
	{
#ifdef RICH_MPI
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		std::vector<size_t> bad_indeces;
		std::vector<size_t> neigh;
		size_t Nrefine = ToRefine.first.size();
		std::vector<Vector3D> res;
		res.reserve(Nrefine);
		// Do we have a prefred direction?
		if (!ToRefine.second.empty())
		{
			for (size_t i = 0; i < Nrefine; ++i)
			{
				double R = tess.GetWidth(ToRefine.first[i]);
				Vector3D newpoint = tess.GetMeshPoint(ToRefine.first[i]) + 0.25*R*ToRefine.second[i] / fastabs(ToRefine.second[i]);
				res.push_back(newpoint);
			}
		}
		else
		{
			// Find direction by opposite to closest point
			for (size_t i = 0; i < Nrefine; ++i)
			{
				// Find closest neighbor
				Vector3D point = tess.GetMeshPoint(ToRefine.first[i]);
				tess.GetNeighbors(ToRefine.first[i], neigh);
				double min_dist2 = ScalarProd(point - tess.GetMeshPoint(neigh[0]), point - tess.GetMeshPoint(neigh[0]));
				size_t index = 0;
				for (size_t j = 1; j < neigh.size(); ++j)
				{
					double temp = ScalarProd(point - tess.GetMeshPoint(neigh[j]), point - tess.GetMeshPoint(neigh[j]));
					if (temp < min_dist2)
					{
						index = j;
						min_dist2 = temp;
					}
				}
				// Make new point
				Vector3D n_split = (point - tess.GetMeshPoint(neigh[index]));
				n_split = normalize(n_split);
				res.push_back(point + 0.25*tess.GetWidth(ToRefine.first[i])*n_split);
			}
		}
		// Make sure point is inside domian
#ifndef RICH_MPI
		std::pair<Vector3D, Vector3D> bb = tess.GetBoxCoordinates();
#endif
		for (size_t i = 0; i < Nrefine; ++i)
		{
#ifdef RICH_MPI
			if (!PointInPoly(tproc, res[i], rank))
				bad_indeces.push_back(i);
#else
			if (res[i].x > bb.second.x || res[i].x<bb.first.x || res[i].y>bb.second.y || res[i].y<bb.first.y
				|| res[i].z>bb.second.z || res[i].z < bb.first.z)
				bad_indeces.push_back(i);
#endif
		}
		if (!bad_indeces.empty())
		{
			RemoveVector(res, bad_indeces);
			RemoveVector(ToRefine.first, bad_indeces);
		}
		return res;
	}

#ifdef RICH_MPI
	void SendRecvMPIFullRemove(Tessellation3D const& tess, vector<size_t> const& toremove, vector<vector<size_t> >
		&nghost_index, vector<vector<vector<size_t> > > &duplicated_index, vector<vector<vector<Vector3D> > > &planes,
		vector<vector<vector<double> > > &planes_d)
	{
		vector<size_t> temp;
		size_t Nprocs = tess.GetDuplicatedProcs().size();
		nghost_index.clear();
		nghost_index.resize(Nprocs);
		duplicated_index.clear();
		duplicated_index.resize(Nprocs);
		planes_d.clear();
		planes_d.resize(Nprocs);
		planes.clear();
		planes.resize(Nprocs);
		vector<vector<size_t> > duplicated_points = tess.GetDuplicatedPoints();
		vector<vector<size_t> > ghost_points = tess.GetGhostIndeces();
		vector<vector<size_t> > sort_indeces(Nprocs), sort_indecesg(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(duplicated_points[i], sort_indeces[i]);
			sort(duplicated_points[i].begin(), duplicated_points[i].end());
			sort_index(ghost_points[i], sort_indecesg[i]);
			sort(ghost_points[i].begin(), ghost_points[i].end());
		}

		size_t nremove = toremove.size();
		size_t Norg = tess.GetPointNo();
		vector<r3d_plane> r_planes;
		for (size_t i = 0; i < nremove; ++i)
		{
			vector<vector<size_t> > ghosts(Nprocs);
			tess.GetNeighbors(toremove[i], temp);
			size_t Nneigh = temp.size();
			bool added = false;
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (temp[j] >= Norg)
				{
					for (size_t k = 0; k < Nprocs; ++k)
					{
						vector<size_t>::const_iterator it = binary_find(ghost_points[k].begin(),
							ghost_points[k].end(), temp[j]);
						if (it != ghost_points[k].end())
						{
							ghosts[k].push_back(sort_indecesg[k][static_cast<size_t>(it - ghost_points[k].begin())]);
							added = true;
						}
					}
				}
			}
			if (added)
			{
				GetPlanes(r_planes, tess, toremove[i]);
				size_t nplanes = r_planes.size();
				vector<Vector3D> v_plane(nplanes);
				vector<double> d_plane(nplanes);
				for (size_t j = 0; j < nplanes; ++j)
				{
					d_plane[j] = r_planes[j].d;
					v_plane[j].x = r_planes[j].n.xyz[0];
					v_plane[j].y = r_planes[j].n.xyz[1];
					v_plane[j].z = r_planes[j].n.xyz[2];
				}
				for (size_t k = 0; k < Nprocs; ++k)
					if (!ghosts[k].empty())
					{
						duplicated_index[k].push_back(ghosts[k]);
						nghost_index[k].push_back(sort_indeces[k][static_cast<size_t>(binary_find(duplicated_points[k].begin(),
							duplicated_points[k].end(), toremove[i]) - duplicated_points[k].begin())]);
						planes[k].push_back(v_plane);
						planes_d[k].push_back(d_plane);
					}
			}
		}
		// send/recv the data
		nghost_index = MPI_exchange_data(tess.GetDuplicatedProcs(), nghost_index);
		duplicated_index = MPI_exchange_data(tess, duplicated_index);
		planes = MPI_exchange_data(tess.GetDuplicatedProcs(), planes, tess.GetMeshPoint(0));
		planes_d = MPI_exchange_data(tess, planes_d);
		// convert the data
		for (size_t i = 0; i < Nprocs; ++i)
		{
			size_t size = nghost_index[i].size();
			assert(duplicated_index[i].size() == size);
			for (size_t j = 0; j < size; ++j)
			{
				nghost_index[i][j] = tess.GetGhostIndeces()[i].at(nghost_index[i][j]);
				size_t size2 = duplicated_index[i][j].size();
				for (size_t k = 0; k < size2; ++k)
					duplicated_index[i][j][k] = tess.GetDuplicatedPoints()[i].at(duplicated_index[i][j][k]);
			}
		}
	}

	void SendRecvMPIRefine(Tessellation3D const& tess, std::vector<std::vector<size_t> > const& to_send,
		std::vector<size_t> const& refined_points, Tessellation3D const& oldtess,
		std::vector<std::vector<std::vector<size_t> > > &neigh_index, std::vector<std::vector<double> > &planes,
		std::vector < std::vector<size_t> > &n_planes, std::vector<std::vector<size_t> > &changed_byouter)
	{
		// to_send is the list of outer neighbors for each refine point
		// refined_points is the index of the refined points in new tess
		assert(refined_points.size() == to_send.size());
		vector<vector<size_t> > duplicated_points = oldtess.GetDuplicatedPoints();
		size_t Nprocs = duplicated_points.size();
		vector<vector<size_t> > ghost_points = oldtess.GetGhostIndeces();
		vector<vector<size_t> > sort_indeces(Nprocs), sort_indecesg(Nprocs);
		// sort the indeces
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(duplicated_points[i], sort_indeces[i]);
			sort(duplicated_points[i].begin(), duplicated_points[i].end());
			sort_index(ghost_points[i], sort_indecesg[i]);
			sort(ghost_points[i].begin(), ghost_points[i].end());
		}
		// Create send data
		neigh_index.clear();
		neigh_index.resize(Nprocs);
		planes.clear();
		planes.resize(Nprocs);
		n_planes.clear();
		n_planes.resize(Nprocs);
		changed_byouter.clear();
		changed_byouter.resize(Nprocs);
		size_t nsend = to_send.size();
		vector<r3d_plane> r_planes;
		for (size_t i = 0; i < nsend; ++i)
		{
			r_planes.clear();
			// Get the polyhedra of new cell
			GetPlanes(r_planes, tess, refined_points[i]);
			// find indeces of neighbors
			for (size_t k = 0; k < Nprocs; ++k)
			{
				std::vector<size_t> indeces_toadd;
				for (size_t j = 0; j < to_send[i].size(); ++j)
				{
					vector<size_t>::const_iterator it = binary_find(ghost_points[k].begin(),
						ghost_points[k].end(), to_send[i][j]);
					if (it != ghost_points[k].end())
						indeces_toadd.push_back(sort_indecesg[k][static_cast<size_t>(it - ghost_points[k].begin())]);
				}
				if (!indeces_toadd.empty())
				{
					changed_byouter[k].push_back(refined_points[i]);
					neigh_index[k].push_back(indeces_toadd);
					size_t nplanes = r_planes.size();
					for (size_t j = 0; j < nplanes; ++j)
					{
						planes[k].push_back(r_planes[j].d);
						planes[k].push_back(r_planes[j].n.xyz[0]);
						planes[k].push_back(r_planes[j].n.xyz[1]);
						planes[k].push_back(r_planes[j].n.xyz[2]);
					}
					n_planes[k].push_back(nplanes);
				}
			}
		}
		// send/recv the data
		neigh_index = MPI_exchange_data(oldtess, neigh_index);
		planes = MPI_exchange_data(oldtess.GetDuplicatedProcs(), planes);
		n_planes = MPI_exchange_data(oldtess.GetDuplicatedProcs(), n_planes);
		// convert the data
		for (size_t i = 0; i < Nprocs; ++i)
		{
			size_t size = neigh_index[i].size();
			for (size_t j = 0; j < size; ++j)
			{
				size_t size2 = neigh_index[i][j].size();
				for (size_t k = 0; k < size2; ++k)
				{
					neigh_index[i][j][k] = oldtess.GetDuplicatedPoints()[i].at(neigh_index[i][j][k]);
				}
			}
		}
	}
#endif

	void LocalRemove(Tessellation3D const& oldtess, std::vector<size_t> const& ToRemove, AMRExtensiveUpdater3D const& eu,
		std::vector<ComputationalCell3D> const& cells, EquationOfState const& eos, TracerStickerNames const& tsn,
		Tessellation3D const& tess, std::vector<Conserved3D> &extensives,SpatialReconstruction3D &interp)
	{
		std::vector<size_t> neigh,temp2;
		point_vec temp;
		r3d_poly poly, poly2;
		std::vector<std::vector<int> > i_temp;
		size_t NRemove = ToRemove.size();
		size_t Norg = oldtess.GetPointNo();
		for (size_t i = 0; i < NRemove; ++i)
		{
			oldtess.GetNeighbors(ToRemove[i], neigh);
			size_t Nneigh = neigh.size();
			// Get old poly
			if (GetPoly(oldtess, ToRemove[i], poly, temp, temp2, i_temp))
			{
				for (size_t j = 0; j < Nneigh; ++j)
				{
					if (neigh[j] >= Norg)
						continue;
					// copy poly
					poly2.nverts = poly.nverts;
					for (int k = 0; k < poly2.nverts; ++k)
					{
						for (size_t l = 0; l < 3; ++l)
						{
							poly2.verts[k].pos.xyz[l] = poly.verts[k].pos.xyz[l];
							poly2.verts[k].pnbrs[l] = poly.verts[k].pnbrs[l];
						}
					}

					size_t index_remove = static_cast<size_t>(std::lower_bound(ToRemove.begin(), ToRemove.end(), neigh[j])
						- ToRemove.begin());
					std::pair<bool, std::array<double,4> > dv = PolyhedraIntersection(tess, neigh[j] - index_remove, poly2);
#ifdef RICH_DEBUG
					try
					{
#endif
						if(dv.first)
							extensives[neigh[j] - index_remove] += eu.ConvertPrimitveToExtensive3D(cells[ToRemove[i]], eos, dv.second[0], tsn, interp.GetSlopes()[ToRemove[i]],
								oldtess.GetCellCM(ToRemove[i]), Vector3D(dv.second[1], dv.second[2], dv.second[3]));
#ifdef RICH_DEBUG
				}
				catch (UniversalError &eo)
				{
					eo.AddEntry("Error in LocalRemove", 0);
					eo.AddEntry("Volume", dv.second[0]);
					eo.AddEntry("Current remove", ToRemove[i]);
					eo.AddEntry("Current remove ID", cells[ToRemove[i]].ID);
					eo.AddEntry("New mass", extensives[neigh[j] - index_remove].mass);
					eo.AddEntry("Remove index", i);
					eo.AddEntry("neigh index",j);
					eo.AddEntry("neigh", neigh[j]);
					eo.AddEntry("index_remove", index_remove);
					throw eo;
				}
#endif

				}
			}
		}
	}

#ifdef RICH_MPI
	void MPIRemove(Tessellation3D const& oldtess, Tessellation3D const& tess, std::vector<size_t> const& ToRemove,
		AMRExtensiveUpdater3D const& eu, EquationOfState const& eos, TracerStickerNames const& tsn,
		std::vector<ComputationalCell3D> const& cells, std::vector<Conserved3D> &extensives,SpatialReconstruction3D &interp)
	{
		vector<vector<size_t> > nghost_index;
		vector<vector<vector<size_t> > > duplicate_index;
		vector<vector < vector<Vector3D> > > planes_v;
		vector<vector < vector<double> > > planes_d;
		vector<r3d_plane> planes;
		r3d_poly poly, poly2;
		std::vector<size_t>  temp2;
		point_vec temp;
		std::vector<std::vector<int> > i_temp;
		SendRecvMPIFullRemove(oldtess, ToRemove, nghost_index, duplicate_index, planes_v, planes_d);
		size_t Nproc = nghost_index.size();
		for (size_t i = 0; i < Nproc; ++i)
		{
			size_t NremoveMPI = duplicate_index[i].size();
			for (size_t j = 0; j < NremoveMPI; ++j)
			{
				// build planes for outer point
				size_t Nplane = planes_v[i].at(j).size();
				planes.resize(Nplane);
				for (size_t k = 0; k < Nplane; ++k)
				{
					planes[k].d = planes_d[i][j].at(k);
					planes[k].n.xyz[0] = planes_v[i][j].at(k).x;
					planes[k].n.xyz[1] = planes_v[i][j][k].y;
					planes[k].n.xyz[2] = planes_v[i][j][k].z;
				}
				size_t NChangeLocal = duplicate_index[i][j].size();
				for (size_t k = 0; k < NChangeLocal; ++k)
				{
					size_t index_remove = static_cast<size_t>(std::lower_bound(ToRemove.begin(), ToRemove.end(),
						duplicate_index[i][j].at(k)) - ToRemove.begin());
					if (GetPoly(tess, duplicate_index[i][j][k] - index_remove, poly, temp, temp2, i_temp))
					{
						// copy poly
						poly2.nverts = poly.nverts;
						for (int kk = 0; kk < poly2.nverts; ++kk)
						{
							for (size_t l = 0; l < 3; ++l)
							{
								poly2.verts[kk].pos.xyz[l] = poly.verts[kk].pos.xyz[l];
								poly2.verts[kk].pnbrs[l] = poly.verts[kk].pnbrs[l];
							}
						}
						std::pair<bool, std::array<double,4> > dv = PolyhedraIntersection(oldtess, 0, poly2, &planes);
#ifdef RICH_DEBUG
						try
						{
#endif
						extensives[duplicate_index[i][j][k] - index_remove] += eu.ConvertPrimitveToExtensive3D(
							cells[nghost_index[i][j]], eos, dv.second[0], tsn,interp.GetSlopes()[nghost_index[i][j]],oldtess.GetCellCM(nghost_index[i][j]),
							Vector3D(dv.second[1], dv.second[2], dv.second[3]));
#ifdef RICH_DEBUG
					}
					catch (UniversalError &eo)
					{
						eo.AddEntry("Error in MPIRemove", 0);
						eo.AddEntry("Volume", dv.second[0]);
						eo.AddEntry("Current remove", nghost_index[i][j]);
						eo.AddEntry("Current remove ID", cells[nghost_index[i][j]].ID);
						eo.AddEntry("New mass", extensives[duplicate_index[i][j][k] - index_remove].mass);
						eo.AddEntry("Remove index i", i);
						eo.AddEntry("Remove index j", j);
						eo.AddEntry("Remove index k", k);
						eo.AddEntry("neigh index", j);
						eo.AddEntry("duplicate_index[i][j][k] - index_remove", duplicate_index[i][j][k] - index_remove);
						eo.AddEntry("index_remove", index_remove);
						throw eo;
					}
#endif

					}
				}
			}
		}
	}
#endif

	void LocalRefine(Tessellation3D const& oldtess, Tessellation3D const& tess, std::vector<size_t> const& ToRefine,
		std::vector<ComputationalCell3D> const& cells, EquationOfState const& eos, TracerStickerNames const& tsn,
		AMRExtensiveUpdater3D const&eu, std::vector<Conserved3D> &extensives,SpatialReconstruction3D &interp)
	{
		size_t Nrefine = ToRefine.size();
		std::vector<size_t> neigh, temp2;
		point_vec temp;
		std::vector<std::vector<int> > i_temp;
		r3d_poly poly, poly2;
		boost::container::flat_set<size_t> checked;
		std::stack<size_t> tocheck;
		size_t Norg = tess.GetPointNo() - ToRefine.size();
		size_t Norg2 = oldtess.GetPointNo();
		for (size_t i = 0; i < Nrefine; ++i)
		{
			checked.clear();
			// Get new cell poly
			if (GetPoly(tess, Norg + i, poly, temp, temp2, i_temp))
			{
				// Get neigh to check
				oldtess.GetNeighbors(ToRefine[i], neigh);
				size_t Nneigh = neigh.size();
				for (size_t j = 0; j < Nneigh; ++j)
					tocheck.push(neigh[j]);
				tocheck.push(ToRefine[i]);
				while (!tocheck.empty())
				{
					size_t cur_check = tocheck.top();
					tocheck.pop();
					// did we check this cell yet?
					if (checked.count(cur_check) == 1)
						continue;
					else // keep track of visted cells
						checked.insert(cur_check);
					if (cur_check >= Norg2)
						continue;
					// copy poly
					poly2.nverts = poly.nverts;
					for (int k = 0; k < poly2.nverts; ++k)
					{
						for (size_t l = 0; l < 3; ++l)
						{
							poly2.verts[k].pos.xyz[l] = poly.verts[k].pos.xyz[l];
							poly2.verts[k].pnbrs[l] = poly.verts[k].pnbrs[l];
						}
					}
					// Check intersectrion
					std::pair<bool, std::array<double, 4> > dv = PolyhedraIntersection(oldtess, cur_check, poly2);
					if (dv.first)
					{
						// Remove extensive from neigh cell and add to new cell
#ifdef RICH_DEBUG
						try
						{
#endif
							Conserved3D toadd = eu.ConvertPrimitveToExtensive3D(cells[cur_check], eos, dv.second[0], tsn, interp.GetSlopes()[cur_check],
								oldtess.GetCellCM(cur_check), Vector3D(dv.second[1], dv.second[2], dv.second[3]));
							extensives[cur_check] -= toadd;
							//extensives[Norg2 + i].tracers.resize(toadd.tracers.size());
							extensives[Norg2 + i] += toadd;
							oldtess.GetNeighbors(cur_check, neigh);
							Nneigh = neigh.size();
							for (size_t j = 0; j < Nneigh; ++j)
								tocheck.push(neigh[j]);
#ifdef RICH_DEBUG
						}
						catch (UniversalError &eo)
						{
							eo.AddEntry("Error in LocalRefine", 0);
							eo.AddEntry("Volume", dv.second[0]);
							eo.AddEntry("Current check", cur_check);
							eo.AddEntry("Current check ID",cells[cur_check].ID);
							eo.AddEntry("Old mass", extensives[cur_check].mass);
							eo.AddEntry("Old density", cells[cur_check].density);
							eo.AddEntry("New mass", extensives[Norg + i].mass);
							eo.AddEntry("Refine index", i);
							eo.AddEntry("Norg", Norg2);
							throw eo;
						}
#endif
					}
				}
			}
			else
			{
				extensives[Norg2 + i] = eu.ConvertPrimitveToExtensive3D(cells[ToRefine[i]], eos, tess.GetVolume(Norg + i), tsn, interp.GetSlopes()[0], 
					Vector3D(), Vector3D());
				extensives[ToRefine[i]] -= extensives[Norg2 + i];
				std::cout << "Warning no good poly localrefine" << std::endl;
			}
		}
	}

#ifdef RICH_MPI
	void MPIRefine(Tessellation3D const& oldtess, Tessellation3D const& tess, std::vector<size_t> const& ToRefine,
		AMRExtensiveUpdater3D const& eu, EquationOfState const& eos, TracerStickerNames const& tsn,
		std::vector<ComputationalCell3D> const& cells, std::vector<Conserved3D> &extensives,SpatialReconstruction3D &interp)
	{
		std::vector<size_t> temp, temp2;
		point_vec ptemp;
		size_t Norg = oldtess.GetPointNo();
		// Get outer refine points
		size_t Nrefine = ToRefine.size();
		std::vector<std::vector<size_t> >  to_send;
		std::vector<size_t> refined_points;
		for (size_t i = 0; i < Nrefine; ++i)
		{
			oldtess.GetNeighbors(ToRefine[i], temp);
			temp2.clear();
			for (size_t j = 0; j < temp.size(); ++j)
			{
				if (temp[j] >= Norg && !oldtess.IsPointOutsideBox(temp[j]))
					temp2.push_back(temp[j]);
			}
			if (!temp2.empty())
			{
				refined_points.push_back(tess.GetPointNo() - Nrefine + i);
				to_send.push_back(temp2);
			}
		}
		// send / recv data
		std::vector<std::vector<size_t> > n_planes;
		std::vector<std::vector<std::vector<size_t> > > neigh_index;
		std::vector < std::vector<double> > planes;
		std::vector<std::vector<size_t> > changed_byouter;
		SendRecvMPIRefine(tess, to_send, refined_points, oldtess, neigh_index, planes, n_planes, changed_byouter);
		// Find the intersections
		r3d_poly poly;
		std::vector<r3d_plane> r_planes;
		std::vector<std::vector<int> > i_temp;
		boost::container::flat_set<size_t> checked;
		std::stack<size_t> tocheck;
		size_t Nprocs = neigh_index.size();
		std::vector<std::vector<Conserved3D> > extensive_tosend(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
		{
			extensive_tosend[i].resize(neigh_index[i].size());
			size_t counter = 0;
			for (size_t j = 0; j < neigh_index[i].size(); ++j)
			{
				checked.clear();
				// Create planes for intersection
				r_planes.resize(n_planes[i][j]);
				for (size_t k = 0; k < n_planes[i][j]; ++k)
				{
					r_planes[k].d = planes[i][counter];
					r_planes[k].n.xyz[0] = planes[i][counter + 1];
					r_planes[k].n.xyz[1] = planes[i][counter + 2];
					r_planes[k].n.xyz[2] = planes[i].at(counter + 3);
					counter += 4;
				}
				// Check for intersections
				for (size_t k = 0; k < neigh_index[i][j].size(); ++k)
					tocheck.push(neigh_index[i][j][k]);
				while (!tocheck.empty())
				{
					size_t cur_check = tocheck.top();
					tocheck.pop();
					// did we check this cell yet?
					if (checked.count(cur_check) == 1)
						continue;
					else // keep track of visted cells
						checked.insert(cur_check);
					if (cur_check >= Norg)
						continue;
					if (GetPoly(oldtess, cur_check, poly, ptemp, temp2, i_temp))
					{
						std::pair<bool, std::array<double,4> > dv = PolyhedraIntersection(oldtess, cur_check, poly, &r_planes);
						if (dv.first)
						{
							// add and remove the extensive
#ifdef RICH_DEBUG
							try
							{
#endif
							Conserved3D toadd = eu.ConvertPrimitveToExtensive3D(cells[cur_check], eos, dv.second[0], tsn, interp.GetSlopes()[cur_check],
								oldtess.GetCellCM(cur_check), Vector3D(dv.second[1], dv.second[2], dv.second[3]));
							extensives[cur_check] -= toadd;
							extensive_tosend[i][j] += toadd;
							oldtess.GetNeighbors(cur_check, temp);
							size_t Nneigh = temp.size();
							for (size_t k = 0; k < Nneigh; ++k)
								tocheck.push(temp[k]);
#ifdef RICH_DEBUG
						}
							catch (UniversalError &eo)
							{
								eo.AddEntry("Error in MPIRefine", 0);
								eo.AddEntry("Volume", dv.second[0]);
								eo.AddEntry("Current check", cur_check);
								eo.AddEntry("Current check ID", cells[cur_check].ID);
								eo.AddEntry("Old mass", extensives[cur_check].mass);
								eo.AddEntry("Old density", cells[cur_check].density);
								eo.AddEntry("Refine index i", i);
								eo.AddEntry("Refine index j", j);
								throw eo;
							}
#endif
						}
					}
					else
					{
						int rank = 0;
						MPI_Comm_rank(MPI_COMM_WORLD, &rank);
						std::cout << "warning bad poly in MPIRefine in rank " <<rank<< std::endl;
					}
				}
			}
		}
		extensive_tosend = MPI_exchange_data(oldtess.GetDuplicatedProcs(), extensive_tosend, extensives.at(0));
		size_t Nremove = oldtess.GetPointNo() + Nrefine - tess.GetPointNo();
		for (size_t i = 0; i < Nprocs; ++i)
		{
			for (size_t j = 0; j < extensive_tosend[i].size(); ++j)
				extensives[changed_byouter[i][j] + Nremove] += extensive_tosend[i][j];
		}
	}
#endif

}

AMRCellUpdater3D::~AMRCellUpdater3D(void) {}

AMRExtensiveUpdater3D::~AMRExtensiveUpdater3D(void) {}

Conserved3D SimpleAMRExtensiveUpdater3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
	double volume, TracerStickerNames const& tracerstickernames, Slope3D const& slope, Vector3D const& CMold, Vector3D const& CMnew) const
{
	Conserved3D res;
	Vector3D diff(CMnew);
	diff -= CMold;
	
	ComputationalCell3D cell_temp(cell);
	ComputationalCellAddMult(cell_temp, slope.xderivative, diff.x);
	ComputationalCellAddMult(cell_temp, slope.yderivative, diff.y);
	ComputationalCellAddMult(cell_temp, slope.zderivative, diff.z);
	cell_temp.internal_energy = eos.dp2e(cell_temp.density, cell_temp.pressure, cell_temp.tracers, tracerstickernames.tracer_names);
	const double mass = volume* cell_temp.density;
	res.mass = mass;
	res.internal_energy = cell_temp.internal_energy*mass;
	res.energy = res.internal_energy + 0.5*mass*ScalarProd(cell_temp.velocity, cell_temp.velocity);
	res.momentum = mass* cell_temp.velocity;
	size_t N = cell_temp.tracers.size();
	//res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell_temp.tracers[i] * mass;
	return res;
}

Conserved3D SimpleAMRExtensiveUpdaterSR3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
	double volume, TracerStickerNames const& tracerstickernames, Slope3D const& slope, Vector3D const& CMold, Vector3D const& CMnew) const
{
	Conserved3D res;
	Vector3D diff(CMnew);
	diff -= CMold;
	ComputationalCell3D cell_temp(cell);
	ComputationalCellAddMult(cell_temp, slope.xderivative, diff.x);
	ComputationalCellAddMult(cell_temp, slope.yderivative, diff.y);
	ComputationalCellAddMult(cell_temp, slope.zderivative, diff.z);
	cell_temp.internal_energy = eos.dp2e(cell_temp.density, cell_temp.pressure, cell_temp.tracers, tracerstickernames.tracer_names);
	double gamma = 1.0 / std::sqrt(1 - ScalarProd(cell_temp.velocity, cell_temp.velocity));
	const double mass = volume * cell_temp.density*gamma;
	res.mass = mass;
	size_t N = cell_temp.tracers.size();
	//res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell_temp.tracers[i] * mass;
	double enthalpy = cell_temp.internal_energy;
	if (fastabs(cell_temp.velocity) < 1e-5)
		res.energy = (gamma*enthalpy + 0.5*ScalarProd(cell_temp.velocity, cell_temp.velocity))* mass - cell_temp.pressure*volume;
	else
		res.energy = (gamma*enthalpy + (gamma - 1))* mass - cell_temp.pressure*volume;
	res.momentum = mass * cell_temp.velocity *gamma*(enthalpy + 1);
	return res;
}

SimpleAMRCellUpdater3D::SimpleAMRCellUpdater3D(vector<string> toskip) :toskip_(toskip) {}

ComputationalCell3D SimpleAMRCellUpdater3D::ConvertExtensiveToPrimitve3D(const Conserved3D& extensive, const EquationOfState& eos,
	double volume, ComputationalCell3D const& old_cell, TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < toskip_.size(); ++i)
	{
		if (*safe_retrieve(old_cell.stickers.begin(), tracerstickernames.sticker_names.begin(),
			tracerstickernames.sticker_names.end(), toskip_[i]))
			return old_cell;
	}

	ComputationalCell3D res;
	const double vol_inv = 1.0 / volume;
	res.density = extensive.mass*vol_inv;
	res.velocity = extensive.momentum / extensive.mass;
	res.ID  = old_cell.ID;
	try
	{
		res.pressure = eos.de2p(res.density, extensive.internal_energy / extensive.mass);
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Density", res.density);
		eo.AddEntry("Vx", res.velocity.x);
		eo.AddEntry("Vy", res.velocity.y);
		eo.AddEntry("Vz", res.velocity.z);
		eo.AddEntry("ID", static_cast<double>(res.ID));
		eo.AddEntry("internal energy", extensive.internal_energy / extensive.mass);
		eo.AddEntry("Volume", 1.0 / vol_inv);
		throw eo;
	}
	res.internal_energy = extensive.internal_energy / extensive.mass;
	size_t N = extensive.tracers.size();
//	res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = extensive.tracers[i] / extensive.mass;
	res.stickers = old_cell.stickers;
	return res;
}

SimpleAMRCellUpdaterSR3D::SimpleAMRCellUpdaterSR3D(double G, vector<string> toskip) : G_(G), toskip_(toskip) {}

ComputationalCell3D SimpleAMRCellUpdaterSR3D::ConvertExtensiveToPrimitve3D(const Conserved3D& extensive, const EquationOfState& /*eos*/,
	double volume, ComputationalCell3D const& old_cell, TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < toskip_.size(); ++i)
		if (*safe_retrieve(old_cell.stickers.begin(), tracerstickernames.sticker_names.begin(),
			tracerstickernames.sticker_names.end(), toskip_[i]))
			return old_cell;
		//if (safe_retrieve(old_cell.stickers, tracerstickernames.sticker_names, toskip_[i]))
			//return old_cell;

	double v = GetVelocity(extensive, G_);
	volume = 1.0 / volume;
	ComputationalCell3D res;
	if (res.density < 0)
		throw UniversalError("Negative density in SimpleAMRCellUpdaterSR");
	res.velocity = (fastabs(extensive.momentum)*1e8 < extensive.mass) ? extensive.momentum / extensive.mass : v * extensive.momentum / abs(extensive.momentum);
	double gamma_1 = std::sqrt(1 - ScalarProd(res.velocity, res.velocity));
	res.density = extensive.mass *gamma_1*volume;
	res.stickers = old_cell.stickers;
	//	size_t N = extensive.tracers.size();
//	res.tracers.resize(N);
	for (size_t i = 0; i < extensive.tracers.size(); ++i)
		res.tracers[i] = extensive.tracers[i] / extensive.mass;
	if (fastabs(res.velocity) < 1e-5)
		res.pressure = (G_ - 1)*((extensive.energy - ScalarProd(extensive.momentum, res.velocity))*volume
			+ (0.5*ScalarProd(res.velocity, res.velocity))*res.density);
	else
		res.pressure = (G_ - 1)*(extensive.energy*volume - ScalarProd(extensive.momentum, res.velocity)*volume
			+ (1.0 / gamma_1 - 1)*res.density);
	return res;
}

SimpleAMRExtensiveUpdater3D::SimpleAMRExtensiveUpdater3D(void) {}

CellsToRemove3D::~CellsToRemove3D(void) {}

CellsToRefine3D::~CellsToRefine3D(void) {}

AMR3D::AMR3D(EquationOfState const& eos, 
	     CellsToRefine3D const& refine,
	     CellsToRemove3D const& remove,
	     SpatialReconstruction3D &interp,
	     AMRCellUpdater3D* cu,
	     AMRExtensiveUpdater3D* eu):
  eos_(eos), 
  refine_(refine), 
  remove_(remove),
  scu_(SimpleAMRCellUpdater3D()),
  seu_(SimpleAMRExtensiveUpdater3D()), 
  interp_(interp), 
  cu_(cu), 
  eu_(eu)
{
	if (!cu)
		cu_ = &scu_;
	if (!eu)
		eu_ = &seu_;
}


void AMR3D::operator() (HDSim3D &sim)
{
	Tessellation3D &tess = sim.getTesselation();
	std::vector<ComputationalCell3D> &cells = sim.getCells();
	std::vector<Conserved3D> &extensives = sim.getExtensives();
	EquationOfState const& eos = eos_;
	TracerStickerNames tsn = sim.GetTracerStickerNames();
	double time = sim.GetTime();
	// Get remove list
	std::pair<vector<size_t>, vector<double> > ToRemove = remove_.ToRemove(tess, cells, time, tsn);
	// sort
	vector<size_t> indeces = sort_index(ToRemove.first);
	ToRemove.second = VectorValues(ToRemove.second, indeces);
	ToRemove.first = VectorValues(ToRemove.first, indeces);
	// remove neighboring remove points
	ToRemove = RemoveNeighbors(ToRemove.second, ToRemove.first, tess);
#ifdef RICH_MPI
	ToRemove = RemoveMPINeighbors(ToRemove.second, ToRemove.first, tess);
#endif
	// Get points to refine
	std::pair<vector<size_t>, std::vector<Vector3D> > ToRefine = refine_.ToRefine(tess, cells, time, tsn);
	sort_index(ToRefine.first, indeces);
	sort(ToRefine.first.begin(), ToRefine.first.end());
	if (!ToRefine.second.empty())
		VectorValues(ToRefine.second, indeces);
	// remove neighboring refine/remove
	RemoveRefineNeighborRemove(tess, ToRemove.first, ToRefine.first, ToRefine.second);

	// Do we need to rebuild tess ?
	int ntotal = static_cast<int>(ToRemove.first.size() + ToRefine.first.size());
#ifdef RICH_MPI
	int ntemp = 0;
	MPI_Allreduce(&ntotal, &ntemp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	ntotal = ntemp;
#endif
	if (ntotal == 0)
		return;
	interp_.BuildSlopes(tess, cells, time, tsn);
	// Get new points from refine
	std::vector<Vector3D> new_points = GetNewPoints(tess, ToRefine
#ifdef RICH_MPI
		, sim.getProcTesselation()
#endif
	);
	// Create copy of old tess
	boost::scoped_ptr<Tessellation3D> oldtess(tess.clone());
	// Build new tess
	vector<Vector3D> new_mesh = tess.GetMeshPoints();
	size_t Norg = tess.GetPointNo();
	new_mesh.resize(Norg);
	RemoveVector(new_mesh, ToRemove.first);
	new_mesh.insert(new_mesh.end(), new_points.begin(), new_points.end());
#ifdef RICH_MPI
	tess.Build(new_mesh, sim.getProcTesselation());
#else
	tess.Build(new_mesh);
#endif
	// Fix extensives for refine
	extensives.resize(oldtess->GetPointNo() + ToRefine.first.size());
	LocalRefine(*oldtess, tess, ToRefine.first, cells, eos, tsn, *eu_, extensives,interp_);
#ifdef RICH_MPI
	MPIRefine(*oldtess, tess, ToRefine.first, *eu_, eos, tsn, cells, extensives,interp_);
#endif
	// Remove from extensive the remove cells
	RemoveVector(extensives, ToRemove.first);
	// Add the removed extensive to the neighboring cells
	LocalRemove(*oldtess, ToRemove.first, *eu_, cells, eos_, tsn, tess, extensives,interp_);

#ifdef RICH_MPI
	MPIRemove(*oldtess, tess, ToRemove.first, *eu_, eos, tsn, cells, extensives,interp_);
#endif
	// Recalc cells
	RemoveVector(cells, ToRemove.first);
	size_t NorgNew = tess.GetPointNo();
	cells.resize(NorgNew);
	for (size_t i = 0; i < (Norg - ToRemove.first.size()); ++i)
	{
		try
		{
			cells[i] = cu_->ConvertExtensiveToPrimitve3D(extensives[i], eos, tess.GetVolume(i), cells[i], tsn);
		}
		catch (UniversalError & eo)
		{
			eo.AddEntry("First loop", static_cast<double>(i));
			eo.AddEntry("Norg", static_cast<double>(Norg));
			throw eo;
		}
	}

	// Get index for ID
	size_t Nrefine = ToRefine.first.size();
	size_t Nstart = sim.GetMaxID() + 1;
#ifdef RICH_MPI
	int ws = 0, rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<size_t> nrecv(static_cast<size_t>(ws), 0);
	MPI_Allgather(&Nrefine, 1, MPI_UNSIGNED_LONG_LONG, &nrecv[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
	for (size_t i = 0; i < static_cast<size_t>(rank); ++i)
		Nstart += nrecv[i];
#endif
	// Add new refined cells	
	for (size_t i = 0; i< ToRefine.first.size(); ++i)
	{
		size_t index_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(),
			ToRefine.first[i]) - ToRemove.first.begin());
		try
		{
			cells[(Norg - ToRemove.first.size()) + i] = cu_->ConvertExtensiveToPrimitve3D(extensives[(Norg -
				ToRemove.first.size()) + i], eos, tess.GetVolume((Norg - ToRemove.first.size()) + i),
				cells[ToRefine.first[i] - index_remove], tsn);
			// Add new ID
			cells[(Norg - ToRemove.first.size()) + i].ID = Nstart + i;
		}
		catch (UniversalError & eo)
		{
			eo.AddEntry("Second loop", static_cast<double>(i));
			eo.AddEntry("Norg", static_cast<double>(Norg));
			eo.AddEntry("Nrefine", static_cast<double>(ToRefine.first.size()));
			throw eo;
		}
	}
	// Recalc entropy if needed
	size_t entropy_index = static_cast<size_t>(std::find(tsn.tracer_names.begin(), tsn.tracer_names.end(), std::string("Entropy")) -
		tsn.tracer_names.begin());
	if (entropy_index < tsn.tracer_names.size())
	{
		size_t Nentropy = cells.size();
		for (size_t i = 0; i < Nentropy; ++i)
		{
			cells[i].tracers[entropy_index] = eos.dp2s(cells[i].density, cells[i].pressure, cells[i].tracers, tsn.tracer_names);
			extensives[i].tracers[entropy_index] = cells[i].tracers[entropy_index] * extensives[i].mass;
		}
	}

#ifdef RICH_MPI
	// Update cells
	ComputationalCell3D cdummy;
	MPI_exchange_data(tess, cells, true,&cdummy);
#endif
	// Update Max ID
	size_t & MaxID = sim.GetMaxID();
#ifdef RICH_MPI
	for (size_t i = 0; i < static_cast<size_t>(ws); ++i)
		MaxID += nrecv[i];
#else
	MaxID += Nrefine;
#endif
}

AMR3D::~AMR3D(void) {}
