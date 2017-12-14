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
#ifdef RICH_MPI
	vector<size_t> GetMPIRefineSend(Tessellation3D const& tess, size_t Nsplit, size_t Ntotal)
	{
		vector<size_t> neigh;
		vector<size_t> to_send;
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < Nsplit; ++i)
		{
			tess.GetNeighbors(Ntotal + i, neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] > Norg && neigh[j] < Ntotal && !tess.IsPointOutsideBox(neigh[j]))
				{
					to_send.push_back(Ntotal + i);
					break;
				}
			}
		}
		return to_send;
	}

	vector<vector<size_t> > SendMPIRefine(Tessellation3D const& tess, vector<size_t> const& tosend, vector<vector<Vector3D> > &
		recv_points, vector<vector<size_t> >&splitted_points, vector<vector<vector<int> > > &recv_neigh, size_t Ntotal0,
		vector<size_t> const& ToRefine)
	{
		size_t send_size = tosend.size();
		splitted_points.clear();
		vector<size_t> neigh;
		vector<vector<size_t> > Nghost = tess.GetGhostIndeces();
		vector<vector<size_t> > duplicated_points = tess.GetDuplicatedPoints();
		size_t Nprocs = Nghost.size();
		vector<vector<size_t> > sort_indeces(Nprocs), sort_indeces_duplicated(Nprocs), sent_points(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(Nghost[i], sort_indeces[i]);
			std::sort(Nghost[i].begin(), Nghost[i].end());
			sort_index(duplicated_points[i], sort_indeces_duplicated[i]);
			std::sort(duplicated_points[i].begin(), duplicated_points[i].end());
		}
		size_t Norg = tess.GetPointNo();
		// Create send data
		vector<vector<Vector3D> > sendpoints(Nprocs);
		vector<vector<vector<int> > > sendNghost(Nprocs);
		splitted_points.resize(Nprocs);
		for (size_t i = 0; i < send_size; ++i)
		{
			tess.GetNeighbors(tosend[i], neigh);
			size_t Nneigh = neigh.size();
			bool good = false;
			for (size_t k = 0; k < Nprocs; ++k)
			{
				vector<int> int_temp;
				for (size_t j = 0; j < Nneigh; ++j)
				{
					if (neigh[j] < Norg || tess.IsPointOutsideBox(neigh[j]))
						continue;
					vector<size_t>::const_iterator it = binary_find(Nghost[k].begin(), Nghost[k].end(), neigh[j]);
					if (it != Nghost[k].end())
					{
						good = true;
						int_temp.push_back(static_cast<int>(sort_indeces[k][static_cast<size_t>(it - Nghost[k].begin())]));
					}
				}
				if (!int_temp.empty())
				{
					vector<size_t>::const_iterator it = binary_find(duplicated_points[k].begin(), duplicated_points[k].end(),
						ToRefine[tosend[i] - Ntotal0]);
					splitted_points[k].push_back(sort_indeces_duplicated[k][static_cast<size_t>(it - duplicated_points[k].begin())]);
					sendpoints[k].push_back(tess.GetMeshPoint(tosend[i]));
					sent_points[k].push_back(tosend[i] - Ntotal0 + Norg);
					sendNghost[k].push_back(int_temp);
				}
			}
			assert(good);
		}
		// Communicate data
		splitted_points = MPI_exchange_data(tess.GetDuplicatedProcs(), splitted_points);
		for (size_t i = 0; i < Nprocs; ++i)
			for (size_t j = 0; j < splitted_points.at(i).size(); ++j)
				splitted_points[i][j] = tess.GetGhostIndeces()[i][splitted_points[i][j]];
		recv_points = MPI_exchange_data(tess.GetDuplicatedProcs(), sendpoints, tess.GetMeshPoint(0));
		recv_neigh = MPI_exchange_data(tess, sendNghost);
		for (size_t i = 0; i < recv_neigh.size(); ++i)
			for (size_t j = 0; j < recv_neigh[i].size(); ++j)
				for (size_t k = 0; k < recv_neigh[i][j].size(); ++k)
					recv_neigh[i][j][k] = static_cast<int>(tess.GetDuplicatedPoints()[i][recv_neigh[i][j][k]]);
		return sent_points;
	}
#endif

	void RemoveBadAspectRatio(Tessellation3D const& tess, std::pair<vector<size_t>, vector<Vector3D> > &toremove)
	{
		std::pair<vector<size_t>, vector<Vector3D > > remove_res;
		size_t N = toremove.first.size();
		for (size_t i = 0; i < N; ++i)
		{
			vector<size_t> const& faces = tess.GetCellFaces(toremove.first[i]);
			size_t Nfaces = faces.size();
			bool good = true;
			double R = tess.GetWidth(toremove.first[i]);
			for (size_t j = 0; j < Nfaces; ++j)
			{
				std::pair<size_t, size_t> fneigh = tess.GetFaceNeighbors(faces[j]);
				if (abs(tess.GetMeshPoint(fneigh.first) - tess.GetMeshPoint(fneigh.second)) < 0.2*R)
				{
					good = false;
					break;
				}
			}
			if (good)
			{
				remove_res.first.push_back(toremove.first[i]);
				if (!toremove.second.empty())
					remove_res.second.push_back(toremove.second[i]);
			}
		}
		toremove = remove_res;
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

	void SendRecvMPIRemoveData(Tessellation3D &tess, vector<size_t> const& to_remove, vector<vector<vector<size_t> > >
		&nghost_neigh_index, vector<vector<size_t> > &nghost_remove, vector<vector<vector<size_t> > > & duplicate_neigh_index,
		vector<vector<size_t> > &local_duplicate_remove, vector<vector<size_t> > &all_remove)
	{
		vector<size_t> temp;
		size_t Nprocs = tess.GetDuplicatedProcs().size();
		nghost_neigh_index.resize(Nprocs);
		duplicate_neigh_index.resize(Nprocs);
		all_remove.clear();
		all_remove.resize(Nprocs);
		vector<vector<Vector3D> > to_add_points(Nprocs), to_add_cm(Nprocs);
		nghost_remove.resize(Nprocs);
		vector<vector<size_t> > duplicated_points = tess.GetDuplicatedPoints();
		vector<vector<size_t> > ghost_points = tess.GetGhostIndeces();
		vector<vector<size_t> > sort_indeces(Nprocs), new_send(Nprocs), sort_indecesg(Nprocs);
		local_duplicate_remove.clear();
		local_duplicate_remove.resize(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(duplicated_points[i], sort_indeces[i]);
			sort(duplicated_points[i].begin(), duplicated_points[i].end());
			sort_index(ghost_points[i], sort_indecesg[i]);
			sort(ghost_points[i].begin(), ghost_points[i].end());
		}

		size_t nremove = to_remove.size();
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < nremove; ++i)
		{
			for (size_t k = 0; k < Nprocs; ++k)
			{
				vector<size_t>::const_iterator it = binary_find(duplicated_points[k].begin(), duplicated_points[k].end(), to_remove[i]);
				if (it != duplicated_points[k].end())
					all_remove[k].push_back(sort_indeces[k][static_cast<size_t>(it - duplicated_points[k].begin())]);
			}
			tess.GetNeighbors(to_remove[i], temp);
			size_t Nneigh = temp.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (temp[j] >= Norg)
				{
					for (size_t k = 0; k < Nprocs; ++k)
					{
						vector<size_t>::const_iterator it = binary_find(duplicated_points[k].begin(),
							duplicated_points[k].end(), to_remove[i]);
						if (it != duplicated_points[k].end())
						{
							vector<size_t> temp2;
							for (size_t z = 0; z < Nneigh; ++z)
							{
								vector<size_t>::const_iterator it2 = binary_find(ghost_points[k].begin(),
									ghost_points[k].end(), temp[z]);
								if (it2 != ghost_points[k].end())
									temp2.push_back(sort_indecesg[k][static_cast<size_t>(it2 - ghost_points[k].begin())]);
							}
							if (temp2.empty())
								continue;
							for (size_t z = 0; z < Nneigh; ++z)
							{
								if (temp[z] < Norg)
								{
									vector<size_t>::const_iterator it3 = binary_find(duplicated_points[k].begin(),
										duplicated_points[k].end(), temp[z]);
									if (it3 == duplicated_points[k].end())
										new_send[k].push_back(temp[z]);
								}
							}
							nghost_remove[k].push_back(sort_indeces[k][static_cast<size_t>(it - duplicated_points[k].begin())]);
							duplicate_neigh_index[k].push_back(temp2);
							local_duplicate_remove[k].push_back(to_remove[i]);
						}
					}
					break;
				}
			}
		}
		// Make Unique the new points
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort(new_send[i].begin(), new_send[i].end());
			new_send[i] = unique(new_send[i]);
		}
		// Add to duplicate points, re-sort them and create send data
		for (size_t i = 0; i < Nprocs; ++i)
		{
			if (new_send[i].empty())
				continue;
			to_add_points[i] = VectorValues(tess.GetMeshPoints(), new_send[i]);
			to_add_cm[i] = VectorValues(tess.GetAllCM(), new_send[i]);
			vector<vector<size_t> > &dp = tess.GetDuplicatedPoints();
			dp[i].insert(dp[i].end(), new_send[i].begin(), new_send[i].end());
		}
		duplicated_points = tess.GetDuplicatedPoints();
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(duplicated_points[i], sort_indeces[i]);
			sort(duplicated_points[i].begin(), duplicated_points[i].end());
		}
		// Create send data of neighboring duplicated points
		for (size_t i = 0; i < nremove; ++i)
		{
			tess.GetNeighbors(to_remove[i], temp);
			size_t Nneigh = temp.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (temp[j] >= Norg)
				{
					for (size_t k = 0; k < Nprocs; ++k)
					{
						vector<size_t>::const_iterator it = binary_find(duplicated_points[k].begin(),
							duplicated_points[k].end(), to_remove[i]);
						if (it != duplicated_points[k].end())
						{
							vector<size_t> temp2;
							for (size_t z = 0; z < Nneigh; ++z)
							{
								vector<size_t>::const_iterator it2 = binary_find(ghost_points[k].begin(),
									ghost_points[k].end(), temp[z]);
								if (it2 != ghost_points[k].end())
									temp2.push_back(sort_indecesg[k][static_cast<size_t>(it2 - ghost_points[k].begin())]);
							}
							if (temp2.empty())
								continue;
							temp2.clear();
							for (size_t z = 0; z < Nneigh; ++z)
							{
								vector<size_t>::const_iterator it3 = binary_find(duplicated_points[k].begin(),
									duplicated_points[k].end(), temp[z]);
								if (it3 != duplicated_points[k].end())
									temp2.push_back(sort_indeces[k][static_cast<size_t>(it3 - duplicated_points[k].begin())]);
							}
							nghost_neigh_index[k].push_back(temp2);
						}
					}
					break;
				}
			}
		}

		// Send/Recv the data
		all_remove = MPI_exchange_data(tess.GetDuplicatedProcs(), all_remove);
		nghost_neigh_index = MPI_exchange_data(tess, nghost_neigh_index);
		duplicate_neigh_index = MPI_exchange_data(tess, duplicate_neigh_index);
		to_add_points = MPI_exchange_data(tess.GetDuplicatedProcs(), to_add_points, tess.GetMeshPoint(0));
		to_add_cm = MPI_exchange_data(tess.GetDuplicatedProcs(), to_add_cm, tess.GetMeshPoint(0));
		nghost_remove = MPI_exchange_data(tess.GetDuplicatedProcs(), nghost_remove);
		// Add the recv points
		for (size_t i = 0; i < Nprocs; ++i)
		{
			for (size_t j = 0; j < to_add_points[i].size(); ++j)
			{
				tess.GetGhostIndeces()[i].push_back(tess.GetMeshPoints().size());
				tess.GetMeshPoints().push_back(to_add_points[i][j]);
				tess.GetAllCM().push_back(to_add_cm[i][j]);
			}
			for (size_t j = 0; j < all_remove[i].size(); ++j)
				all_remove[i][j] = tess.GetGhostIndeces()[i][all_remove[i][j]];
		}
	}
#endif //RICH_MPI

	Vector3D GetNewPoint(Tessellation3D const& tess, vector<size_t> const& neigh, size_t index, std::pair<vector<size_t>,
		vector<Vector3D> > &ToRefine)
	{
		if (!ToRefine.second.empty())
		{
			size_t loc = ToRefine.first[index];
			return (tess.GetMeshPoint(loc) + (3e-5*tess.GetWidth(loc) / abs(ToRefine.second.at(index)))*ToRefine.second.at(index));
		}
		size_t Nneigh = neigh.size();
		assert(Nneigh > 0);
		index = ToRefine.first[index];
		Vector3D const& point = tess.GetMeshPoint(index);
		Vector3D other = tess.GetMeshPoint(neigh[0]);
		double max_dist = ScalarProd(point - other, point - other);
		size_t max_loc = 0;
		for (size_t i = 1; i < Nneigh; ++i)
		{
			other = tess.GetMeshPoint(neigh[i]);
			double temp = ScalarProd(point - other, point - other);
			if (temp > max_dist)
			{
				max_dist = temp;
				max_loc = i;
			}
		}
		double eps = 2e-5;
		return point*(1 - eps) + eps*tess.GetMeshPoint(neigh[max_loc]);
	}

	void BuildLocalVoronoi(Tessellation3D &local, Tessellation3D const& tess, vector<size_t> const& real_points,
		Vector3D const& newpoint, size_t torefine)
	{
		size_t Nreal = real_points.size();
		vector<Vector3D> points, ghost;
		points.push_back(newpoint);
		points.push_back(tess.GetMeshPoint(torefine));
		ghost.reserve(Nreal);
		for (size_t i = 0; i < Nreal; ++i)
			ghost.push_back(tess.GetMeshPoint(real_points[i]));
		vector<size_t> v_duplicate(1, 0);
		local.BuildNoBox(points, ghost, v_duplicate);
	}

#ifdef RICH_MPI
	void BuildLocalVoronoiMPI(Tessellation3D &local, Tessellation3D const& tess, vector<int> const& neigh_points,
		Vector3D const& newpoint, size_t refined_index)
	{
		vector<size_t> temp, temp2;
		vector<Vector3D> points, ghost;
		for (size_t i = 0; i < neigh_points.size(); ++i)
		{
			tess.GetNeighbors(static_cast<size_t>(neigh_points[i]), temp);
			temp2.insert(temp2.end(), temp.begin(), temp.end());
			points.push_back(tess.GetMeshPoint(static_cast<size_t>(neigh_points[i])));
		}
		sort(temp2.begin(), temp2.end());
		temp2 = unique(temp2);
		RemoveVal(temp2, refined_index);
		for (size_t i = 0; i < neigh_points.size(); ++i)
			RemoveVal(temp2, static_cast<size_t>(neigh_points[i]));
		size_t Nghost = temp2.size();
		ghost.reserve(Nghost + 1);
		ghost.push_back(newpoint);
		ghost.push_back(tess.GetMeshPoint(refined_index));
		for (size_t i = 0; i < Nghost; ++i)
			ghost.push_back(tess.GetMeshPoint(temp2[i]));
		local.BuildNoBox(points, ghost, vector<size_t>());
	}
#endif

	void FixVoronoi(Tessellation3D &local, Tessellation3D &tess, vector<size_t> &neigh, size_t torefine,
		double &vol, Vector3D &CM, size_t Ntotal0, size_t index)
	{
		// neigh is sorted
		vector<std::pair<size_t, size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		vector<double> & full_area = tess.GetAllArea();
		vector<Vector3D> &full_face_cm = tess.GetAllFaceCM();
		vector<Vector3D>& full_vertices = tess.GetFacePoints();

		size_t Norg = tess.GetPointNo();
		vector<size_t> faces = tess.GetCellFaces(torefine);
		size_t Nfaces = faces.size();
		// Remove old face reference
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t other = tess.GetFaceNeighbors(faces[i]).first == torefine ? tess.GetFaceNeighbors(faces[i]).second :
				tess.GetFaceNeighbors(faces[i]).first;
			if (other < Norg || (other >= Ntotal0 && other < (Ntotal0 + index)))
				RemoveVal(full_cellfaces[other], faces[i]);
		}
		full_cellfaces[torefine].clear();

		tess.GetAllVolumes()[torefine] = local.GetVolume(1);
		tess.GetAllCM()[torefine] = local.GetCellCM(1);
		vol = local.GetVolume(0);
		CM = local.GetCellCM(0);

		// Add new faces
		vector<size_t> temp;
		std::pair<size_t, size_t> new_face_neigh;
		faces = local.GetCellFaces(0);
		faces.insert(faces.end(), local.GetCellFaces(1).begin(), local.GetCellFaces(1).end());
		sort(faces.begin(), faces.end());
		faces = unique(faces);
		Nfaces = faces.size();
		size_t Nlocal = neigh.size() + 6;
		size_t Nvert = full_vertices.size();
		full_cellfaces.resize(Ntotal0 + index + 1);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (localfaceneigh[faces[i]].first == 0)
				new_face_neigh.first = Ntotal0 + index;
			else
				if (localfaceneigh[faces[i]].first == 1)
					new_face_neigh.first = torefine;
				else
					new_face_neigh.first = neigh[localfaceneigh[faces[i]].first - 6];
			if (localfaceneigh[i].second == 1)
				new_face_neigh.second = torefine;
			else
			{
				if (localfaceneigh[faces[i]].second < Nlocal)
					new_face_neigh.second = neigh[localfaceneigh[faces[i]].second - 6];
				else
				{ // New boundary point
					new_face_neigh.second = tess.GetMeshPoints().size();
					tess.GetMeshPoints().push_back(local.GetMeshPoint(localfaceneigh[faces[i]].second));
					tess.GetAllCM().push_back(local.GetMeshPoint(localfaceneigh[faces[i]].second) + local.GetMeshPoint(0)
						- CM);
				}
			}
			temp = local.GetPointsInFace(faces[i]);
			size_t N = temp.size();
			for (size_t j = 0; j < N; ++j)
				temp[j] += Nvert;
			if (new_face_neigh.second < new_face_neigh.first)
			{
				size_t ttemp = new_face_neigh.second;
				new_face_neigh.second = new_face_neigh.first;
				new_face_neigh.first = ttemp;
				FlipVector(temp);
			}
			full_facepoints.push_back(temp);
			full_faceneigh.push_back(new_face_neigh);
			if (new_face_neigh.first < Norg || (new_face_neigh.first >= Ntotal0 && new_face_neigh.first <= (Ntotal0 + index)))
				full_cellfaces.at(new_face_neigh.first).push_back(full_faceneigh.size() - 1);
			if (new_face_neigh.second < Norg || ((new_face_neigh.second >= Ntotal0) && (new_face_neigh.second <= Ntotal0 + index)))
				full_cellfaces.at(new_face_neigh.second).push_back(full_faceneigh.size() - 1);
			full_area.push_back(local.GetArea(i));
			full_face_cm.push_back(local.FaceCM(i));
		}
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}
#ifdef RICH_MPI
	void FixVoronoiMPI(Tessellation3D &local, Tessellation3D &tess, size_t refinedghost,
		vector<int> const& refined_neigh, vector<size_t> &bad_faces)
	{
		// neigh is sorted
		vector<std::pair<size_t, size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		vector<double> & full_area = tess.GetAllArea();
		vector<Vector3D> &full_face_cm = tess.GetAllFaceCM();
		vector<Vector3D>& full_vertices = tess.GetFacePoints();

		size_t Nrefine_neigh = refined_neigh.size();
		// Remove old face reference
		for (size_t i = 0; i < Nrefine_neigh; ++i)
		{
			vector<size_t> faces = tess.GetCellFaces(static_cast<size_t>(refined_neigh[i]));
			for (size_t j = 0; j < faces.size(); ++j)
			{
				if (tess.GetFaceNeighbors(faces[j]).second == refinedghost)
				{
					RemoveVal(full_cellfaces[static_cast<size_t>(refined_neigh[i])], faces[j]);
					bad_faces.push_back(faces[j]);
					break;
				}
			}
		}

		// Add new faces
		vector<size_t> temp, temp2;
		size_t Nvert = full_vertices.size();
		std::pair<size_t, size_t> new_face_neigh;
		for (size_t i = 0; i < Nrefine_neigh; ++i)
		{
			temp = local.GetCellFaces(i);
			size_t Nfaces = temp.size();
			new_face_neigh.first = static_cast<size_t>(refined_neigh[i]);
			for (size_t j = 0; j < Nfaces; ++j)
			{
				bool good = false;
				if (localfaceneigh[temp[j]].second == (Nrefine_neigh + 4))
				{
					new_face_neigh.second = tess.GetMeshPoints().size() - 1;
					good = true;
				}
				if (localfaceneigh[temp[j]].second == (Nrefine_neigh + 5))
				{
					new_face_neigh.second = refinedghost;
					good = true;
				}
				if (good)
				{
					temp2 = local.GetPointsInFace(temp[j]);
					size_t N = temp2.size();
					for (size_t k = 0; k < N; ++k)
						temp2[k] += Nvert;
					full_facepoints.push_back(temp2);
					full_faceneigh.push_back(new_face_neigh);
					full_cellfaces.at(new_face_neigh.first).push_back(full_faceneigh.size() - 1);
					full_area.push_back(local.GetArea(temp[j]));
					full_face_cm.push_back(local.FaceCM(temp[j]));
				}
			}
		}
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}

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
#endif

	void FixBadIndeces(Tessellation3D &tess, vector<size_t> const& bad_indeces, size_t Nsplit, size_t Ntotal0)
	{
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < Norg; ++i)
		{
			vector<size_t> & faces = tess.GetAllCellFaces()[i];
			size_t nfaces = faces.size();
			for (size_t j = 0; j < nfaces; ++j)
			{
				size_t toremove = static_cast<size_t>(std::lower_bound(bad_indeces.begin(), bad_indeces.end(), faces[j])
					- bad_indeces.begin());
				faces[j] -= toremove;
			}
		}

		size_t Nfaces = tess.GetAllFaceNeighbors().size();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (tess.GetFaceNeighbors(i).first >= Ntotal0)
			{
				if (tess.GetFaceNeighbors(i).first < Ntotal0 + Nsplit)
					tess.GetAllFaceNeighbors()[i].first += Norg - Ntotal0 - Nsplit;
			}
			else
			{
				if (tess.GetFaceNeighbors(i).first >= (Norg - Nsplit))
					tess.GetAllFaceNeighbors()[i].first += Nsplit;
			}
			if (tess.GetFaceNeighbors(i).second >= Ntotal0)
			{
				if (tess.GetFaceNeighbors(i).second < Ntotal0 + Nsplit)
					tess.GetAllFaceNeighbors()[i].second += Norg - Ntotal0 - Nsplit;
			}
			else
			{
				if (tess.GetFaceNeighbors(i).second >= (Norg - Nsplit))
					tess.GetAllFaceNeighbors()[i].second += Nsplit;
			}
			if (tess.GetFaceNeighbors(i).first > tess.GetFaceNeighbors(i).second)
			{
				size_t temp = tess.GetAllFaceNeighbors()[i].second;
				tess.GetAllFaceNeighbors()[i].second = tess.GetAllFaceNeighbors()[i].first;
				tess.GetAllFaceNeighbors()[i].first = temp;
				FlipVector(tess.GetAllPointsInFace()[i]);
			}
		}
		vector<vector<size_t> > & ghost = tess.GetGhostIndeces();
		size_t Nghost = ghost.size();
		for (size_t i = 0; i < Nghost; ++i)
		{
			size_t Nghost2 = ghost[i].size();
			for (size_t j = 0; j < Nghost2; ++j)
				if (ghost[i][j] < (Ntotal0 + Nsplit))
					ghost[i][j] += Nsplit;
		}
	}
}

AMRCellUpdater3D::~AMRCellUpdater3D(void) {}

AMRExtensiveUpdater3D::~AMRExtensiveUpdater3D(void) {}

Conserved3D SimpleAMRExtensiveUpdater3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& /*eos*/,
	double volume, TracerStickerNames const& /*tracerstickernames*/) const
{
	Conserved3D res;
	const double mass = volume*cell.density;
	res.mass = mass;
	res.internal_energy = cell.internal_energy*mass;
	res.energy = res.internal_energy + 0.5*mass*ScalarProd(cell.velocity, cell.velocity);
	res.momentum = mass*cell.velocity;
	size_t N = cell.tracers.size();
	res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell.tracers[i] * mass;
	return res;
}

SimpleAMRCellUpdater3D::SimpleAMRCellUpdater3D(vector<string> toskip) :toskip_(toskip) {}

ComputationalCell3D SimpleAMRCellUpdater3D::ConvertExtensiveToPrimitve3D(const Conserved3D& extensive, const EquationOfState& eos,
	double volume, ComputationalCell3D const& old_cell, TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < toskip_.size(); ++i)
		if (safe_retrieve(old_cell.stickers, tracerstickernames.sticker_names, toskip_[i]))
			return old_cell;
	ComputationalCell3D res;
	const double vol_inv = 1.0 / volume;
	res.density = extensive.mass*vol_inv;
	res.velocity = extensive.momentum / extensive.mass;
	res.pressure = eos.de2p(res.density, extensive.internal_energy / extensive.mass);
	res.internal_energy = extensive.internal_energy / extensive.mass;
	size_t N = extensive.tracers.size();
	res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = extensive.tracers[i] / extensive.mass;
	res.stickers = old_cell.stickers;
	return res;
}

SimpleAMRExtensiveUpdater3D::SimpleAMRExtensiveUpdater3D(void) {}

CellsToRemove3D::~CellsToRemove3D(void) {}

CellsToRefine3D::~CellsToRefine3D(void) {}

AMR3D::AMR3D(EquationOfState const& eos, CellsToRefine3D const& refine, CellsToRemove3D const& remove, LinearGauss3D *slopes, AMRCellUpdater3D* cu,
	AMRExtensiveUpdater3D* eu) :eos_(eos), refine_(refine), remove_(remove), scu_(SimpleAMRCellUpdater3D()),
	seu_(SimpleAMRExtensiveUpdater3D()), cu_(cu), eu_(eu), interp_(slopes)
{
	if (!cu)
		cu_ = &scu_;
	if (!eu)
		eu_ = &seu_;
}


void AMR3D::UpdateCellsRefine(Tessellation3D &tess, vector<ComputationalCell3D> &cells, EquationOfState const& /*eos*/,
	vector<Conserved3D> &extensives, double time,
#ifdef RICH_MPI
	Tessellation3D const& /*proctess*/,
#endif
	TracerStickerNames const& tracerstickernames)const
{
	size_t Norg = tess.GetPointNo();
	size_t Ntotal0 = tess.GetMeshPoints().size();
	extensives.resize(Norg);
	cells.resize(Norg);
	std::pair<vector<size_t>, vector<Vector3D> > ToRefine = refine_.ToRefine(tess, cells, time, tracerstickernames);
	vector<size_t> indeces;
	sort_index(ToRefine.first, indeces);
	sort(ToRefine.first.begin(), ToRefine.first.end());
	if (!ToRefine.second.empty())
		VectorValues(ToRefine.second, indeces);
	RemoveBadAspectRatio(tess, ToRefine);
#ifndef RICH_MPI
	if (ToRefine.first.empty())
		return;
#endif
	size_t Nsplit = ToRefine.first.size();
	Voronoi3D vlocal(tess.GetBoxCoordinates().first, tess.GetBoxCoordinates().second);
	vector<size_t> neigh, bad_faces, all_bad_faces, refined, newboundary_faces;
	double newvol;
	Vector3D newCM;
	vector<Vector3D> newpoints, newCMs;
	vector<double> newvols;
	tess.GetMeshPoints().resize(Ntotal0 + Nsplit);
	tess.GetAllCellFaces().resize(Ntotal0);
	for (size_t i = 0; i < Nsplit; ++i)
	{
		tess.GetNeighbors(ToRefine.first[i], neigh);
		sort(neigh.begin(), neigh.end());
		Vector3D NewPoint = GetNewPoint(tess, neigh, i, ToRefine);
		BuildLocalVoronoi(vlocal, tess, neigh, NewPoint, ToRefine.first[i]);
		bad_faces = tess.GetCellFaces(ToRefine.first[i]);
		/////////////
		double oldv = tess.GetVolume(ToRefine.first[i]);
		double newv = vlocal.GetVolume(0) + vlocal.GetVolume(1);
		if (oldv < 0.999*newv || newv < 0.999*oldv)
		{
#ifdef RICH_MPI
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
			std::cout << "Warning, refine old volume " << oldv << " new volume " << newv << " split cell " << ToRefine.first[i] << " index " << i << "volumes are " <<
				vlocal.GetVolume(0) << " " << vlocal.GetVolume(1) <<
#ifdef RICH_MPI
				" rank " << rank <<
#endif
				std::endl;
		}
		///////////
		FixVoronoi(vlocal, tess, neigh, ToRefine.first[i], newvol, newCM, Ntotal0, i);
		PrimitiveToConserved(cells[ToRefine.first[i]], tess.GetVolume(ToRefine.first[i]), extensives[ToRefine.first[i]]);
		tess.GetMeshPoints()[Ntotal0 + i] = NewPoint;
		newvols.push_back(newvol);
		newCMs.push_back(newCM);
		newpoints.push_back(NewPoint);
		all_bad_faces.insert(all_bad_faces.end(), bad_faces.begin(), bad_faces.end());
	}
	// MPI
#ifdef RICH_MPI
	vector<vector<Vector3D> > new_points;
	vector<vector<vector<int> > > recv_neigh;
	vector<vector<size_t> > splitted_points;
	vector<size_t> points_tosend = GetMPIRefineSend(tess, Nsplit, Ntotal0);
	vector<vector<size_t> > sent_points = SendMPIRefine(tess, points_tosend, new_points, splitted_points, recv_neigh, Ntotal0,
		ToRefine.first);
	for (size_t i = 0; i < splitted_points.size(); ++i)
	{
		for (size_t j = 0; j < splitted_points[i].size(); ++j)
		{
			tess.GetGhostIndeces()[i].push_back(tess.GetMeshPoints().size());
			tess.GetMeshPoints().push_back(new_points[i][j]);
			BuildLocalVoronoiMPI(vlocal, tess, recv_neigh[i][j], new_points[i][j], splitted_points[i][j]);
			FixVoronoiMPI(vlocal, tess, splitted_points[i][j], recv_neigh[i][j], all_bad_faces);
		}
}
#endif
	sort(all_bad_faces.begin(), all_bad_faces.end());
	all_bad_faces = unique(all_bad_faces);
	vector<double> & allvol = tess.GetAllVolumes();
	allvol.insert(allvol.begin() + Norg, newvols.begin(), newvols.end());
	vector<Vector3D> & allCM = tess.GetAllCM();
	allCM.insert(allCM.begin() + Norg, newCMs.begin(), newCMs.end());
	// do new cells 
	extensives.resize(Norg + Nsplit);
	for (size_t i = 0; i < newCMs.size(); ++i)
	{
		cells.push_back(cells[ToRefine.first[i]]);
		cells.back().stickers.assign(cells.back().stickers.size(), false);
		PrimitiveToConserved(cells[ToRefine.first[i]], tess.GetVolume(Norg + i), extensives[Norg + i]);
	}
	// Fix the old data and insert new data
	vector<Vector3D> & allpoints = tess.GetMeshPoints();
	size_t Ntemp = allpoints.size();
	allpoints.insert(allpoints.begin() + Norg, newpoints.begin(), newpoints.end());
	allpoints.erase(allpoints.begin() + Ntotal0 + Nsplit, allpoints.begin() + Ntotal0 + 2 * Nsplit);
	for (size_t i = 0; i < newboundary_faces.size(); ++i)
		tess.GetAllFaceNeighbors()[newboundary_faces[i]].second = Ntemp + i;
	size_t &Norg2 = tess.GetPointNo();
	Norg2 += newpoints.size();
	vector<vector<size_t> > temp_cell_faces(tess.GetAllCellFaces().begin() + Ntotal0, tess.GetAllCellFaces().begin() + Nsplit + Ntotal0);
	std::copy(temp_cell_faces.begin(), temp_cell_faces.end(), tess.GetAllCellFaces().begin() + Norg);
	tess.GetAllCellFaces().resize(Norg2);

	RemoveVector(tess.GetAllFaceNeighbors(), all_bad_faces);
	RemoveVector(tess.GetAllArea(), all_bad_faces);
	RemoveVector(tess.GetAllPointsInFace(), all_bad_faces);
	RemoveVector(tess.GetAllFaceCM(), all_bad_faces);
	// Fix all the indeces
	FixBadIndeces(tess, all_bad_faces, Nsplit, Ntotal0);

	// Deal with mpi
#ifdef RICH_MPI
	// Update duplicatedpoints
	for (size_t i = 0; i < sent_points.size(); ++i)
		for (size_t j = 0; j < sent_points[i].size(); ++j)
			tess.GetDuplicatedPoints()[i].push_back(sent_points[i][j]);
	// Update cells and CM
	MPI_exchange_data(tess, tess.GetAllCM(), true);
	MPI_exchange_data(tess, cells, true);
#endif


#ifdef debug_amr
	CheckCorrect(tess);
#endif
	}

void AMR3D::UpdateCellsRemove2(Tessellation3D &tess, vector<ComputationalCell3D> &cells, vector<Conserved3D> &extensives,
	EquationOfState const& eos, double time, TracerStickerNames const& tracerstickernames
#ifdef RICH_MPI
	, Tessellation3D const& proctess
#endif
)const
{
	std::pair<vector<size_t>, vector<double> > ToRemove = remove_.ToRemove(tess, cells, time, tracerstickernames);
	vector<size_t> indeces = sort_index(ToRemove.first);
	ToRemove.second = VectorValues(ToRemove.second, indeces);
	ToRemove.first = VectorValues(ToRemove.first, indeces);
	ToRemove = RemoveNeighbors(ToRemove.second, ToRemove.first, tess);
#ifdef RICH_MPI
	int i_toremove = static_cast<int>(ToRemove.first.size());
	int n_recv = 0;
	MPI_Allreduce(&i_toremove, &n_recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	if (n_recv == 0)
		return;
	ToRemove = RemoveMPINeighbors(ToRemove.second, ToRemove.first, tess);
#else
	if (ToRemove.first.size() == 0)
		return;
#endif
	size_t NRemove = ToRemove.first.size();
	boost::scoped_ptr<Tessellation3D> oldtess(tess.clone());
	vector<Vector3D> new_mesh = tess.GetMeshPoints();
	size_t Norg = tess.GetPointNo();
	new_mesh.resize(Norg);
	RemoveVector(new_mesh, ToRemove.first);
#ifdef RICH_MPI
	tess.Build(new_mesh, proctess);
#else
	tess.Build(new_mesh);
#endif
	RemoveVector(extensives, ToRemove.first);
	vector<size_t> neigh, temp, temp2, changed_cells, changed_cells_old;
	vector<vector<int> > i_temp;
	vector<Vector3D> vtemp;
	r3d_poly poly, poly2;
	// deal with local remove
	for (size_t i = 0; i < NRemove; ++i)
	{
		oldtess->GetNeighbors(ToRemove.first[i], neigh);
		size_t Nneigh = neigh.size();
		GetPoly(*oldtess, ToRemove.first[i], poly, temp, temp2, i_temp);
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

			size_t index_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(), neigh[j])
				- ToRemove.first.begin());
			std::pair<bool, double> dv = PolyhedraIntersection(tess, neigh[j] - index_remove, poly2);
			extensives[neigh[j] - index_remove] += eu_->ConvertPrimitveToExtensive3D(cells[ToRemove.first[i]], eos, dv.second,
				tracerstickernames);
			changed_cells.push_back(neigh[j] - index_remove);
			changed_cells_old.push_back(neigh[j]);
		}
	}
	// deal with mpi cells
#ifdef RICH_MPI
	vector<vector<size_t> > nghost_index;
	vector<vector<vector<size_t> > > duplicate_index;
	vector<vector < vector<Vector3D> > > planes_v;
	vector<vector < vector<double> > > planes_d;
	vector<r3d_plane> planes;
	SendRecvMPIFullRemove(*oldtess, ToRemove.first, nghost_index, duplicate_index, planes_v, planes_d);
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
				size_t index_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(),
					duplicate_index[i][j].at(k)) - ToRemove.first.begin());
				GetPoly(tess, duplicate_index[i][j][k] - index_remove, poly, temp, temp2, i_temp);
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
				std::pair<bool, double> dv = PolyhedraIntersection(*oldtess, 0, poly2, &planes);
				extensives[duplicate_index[i][j][k] - index_remove] += eu_->ConvertPrimitveToExtensive3D(
					cells[nghost_index[i][j]], eos, dv.second, tracerstickernames);
				changed_cells.push_back(duplicate_index[i][j][k] - index_remove);
				changed_cells_old.push_back(duplicate_index[i][j][k]);
			}
		}
	}
#endif
	std::sort(changed_cells_old.begin(), changed_cells_old.end());
	changed_cells_old = unique(changed_cells_old);
	std::sort(changed_cells.begin(), changed_cells.end());
	changed_cells = unique(changed_cells);
	size_t Nchange = changed_cells.size();
	for (size_t i = 0; i < Nchange; ++i)
		cells[changed_cells_old[i]] = cu_->ConvertExtensiveToPrimitve3D(extensives[changed_cells[i]], eos,
			tess.GetVolume(changed_cells[i]), cells[changed_cells_old[i]], tracerstickernames);
	RemoveVector(cells, ToRemove.first);
#ifdef RICH_MPI
	// Update cells and CM
	MPI_exchange_data(tess, cells, true);
#endif
}


void AMR3D::operator() (HDSim3D &sim)
{
	UpdateCellsRemove2(sim.getTesselation(), sim.getCells(), sim.getExtensives(), eos_, sim.GetTime(),
		sim.GetTracerStickerNames()
#ifdef RICH_MPI
		, sim.getProcTesselation()
#endif
	);
	UpdateCellsRefine(sim.getTesselation(), sim.getCells(), eos_, sim.getExtensives(), sim.GetTime(),
#ifdef RICH_MPI
		sim.getProcTesselation(),
#endif
		sim.GetTracerStickerNames());
}

AMR3D::~AMR3D(void) {}