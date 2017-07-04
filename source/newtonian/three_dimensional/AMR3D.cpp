#include "AMR3D.hpp"
#include "../../3D/GeometryCommon/Voronoi3D.hpp"
#include "../../misc/utils.hpp"
#include <boost/array.hpp>
#include <iostream>

//#define debug_amr 1

namespace
{
#ifdef RICH_MPI
	vector<vector<size_t> > GetSentIndeces(Tessellation3D const& tess, vector<size_t> const& ToRemove,
		vector<vector<size_t> > &RemoveIndex)
	{
		RemoveIndex.clear();
		vector<vector<size_t> > sentpoints = tess.GetDuplicatedPoints();
		vector<vector<size_t> > sort_indeces(sentpoints.size());
		vector<vector<size_t> > res(sentpoints.size());
		RemoveIndex.resize(sentpoints.size());
		size_t Nprocs = static_cast<int>(sentpoints.size());
		// sort vectors for fast search
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(sentpoints[i], sort_indeces[i]);
			sort(sentpoints[i].begin(), sentpoints[i].end());
		}
		// search the vectors
		size_t Nremove = ToRemove.size();
		for (size_t i = 0; i < Nremove; ++i)
		{
			for (size_t j = 0; j < Nprocs; ++j)
			{
				vector<size_t>::const_iterator it = binary_find(sentpoints[j].begin(), sentpoints[j].end(),
					ToRemove[i]);
				if (it != sentpoints[j].end())
				{
					res[j].push_back(sort_indeces[j][static_cast<size_t>(it - sentpoints[j].begin())]);
					RemoveIndex[j].push_back(i);
				}
			}
		}
		return res;
	}

	void SendRecvOuterMerits(Tessellation3D const& tess, vector<vector<size_t> > &sent_indeces,
		vector<double> const& merits, vector<vector<size_t> > const& RemoveIndex, vector<vector<size_t> > &recv_indeces,
		vector<vector<double> > &recv_mertis)
	{
		vector<vector<double> > send_merits(sent_indeces.size());
		size_t Nproc = send_merits.size();
		for (size_t i = 0; i < Nproc; ++i)
			send_merits[i] = VectorValues(merits, RemoveIndex[i]);
		recv_indeces = MPI_exchange_data(tess.GetDuplicatedProcs(), sent_indeces);
		recv_mertis = MPI_exchange_data(tess.GetDuplicatedProcs(), send_merits);
	}

	vector<size_t> KeepMPINeighbors(Tessellation3D const& tess, vector<size_t> const& ToRemove,
		vector<double> const& merits, vector<vector<size_t> > &recv_indeces, vector<vector<double> > &recv_mertis)
	{
		vector<vector<size_t> > const& Nghost = tess.GetGhostIndeces();
		size_t Nproc = recv_indeces.size();
		for (size_t i = 0; i < Nproc; ++i)
		{
			recv_indeces[i] = VectorValues(Nghost[i], recv_indeces[i]);
			vector<size_t> temp = sort_index(recv_indeces[i]);
			sort(recv_indeces[i].begin(), recv_indeces[i].end());
			recv_mertis[i] = VectorValues(recv_mertis[i], temp);
		}

		size_t N = tess.GetPointNo();
		size_t Nremove = ToRemove.size();
		vector<size_t> neigh;
		vector<size_t> RemoveFinal;
		for (size_t i = 0; i < Nremove; ++i)
		{
			bool good = true;
			tess.GetNeighbors(ToRemove[i], neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] >= N)
				{
					for (size_t k = 0; k < Nproc; ++k)
					{
						vector<size_t>::const_iterator it = binary_find(recv_indeces[k].begin(), recv_indeces[k].end(),
							neigh[j]);
						if (it != recv_indeces[k].end())
						{
							if (recv_mertis[k][static_cast<size_t>(it - recv_indeces[k].begin())] > merits[i])
							{
								good = false;
								break;
							}
						}
					}
					if (!good)
						break;
				}
			}
			if (good)
				RemoveFinal.push_back(ToRemove[i]);
		}
		return RemoveFinal;
	}

	std::pair<int, size_t> GetOtherCpuInfo(vector<int> const& procnames, vector<vector<size_t> > const&
		duplicated_points, vector<vector<size_t> > const& indeces, size_t point)
	{
		size_t nprocs = procnames.size();
		for (size_t i = 0; i < nprocs; ++i)
		{
			vector<size_t>::const_iterator it = binary_find(duplicated_points[i].begin(), duplicated_points[i].end(),
				point);
			if (it != duplicated_points[i].end())
				return std::pair<int, size_t>(procnames[i], indeces[i][static_cast<size_t>(it - duplicated_points[i].begin())]);
		}
		throw("No matching cpu in amr remove");
	}

	void GetSendData(Tessellation3D &tess, vector<size_t> const& ToRemove,
		vector<Conserved3D> const& extensives, vector<vector<Conserved3D> > &extensive_remove,
		vector<vector<vector<Vector3D> > > &new_point_remove, vector<vector<vector<size_t> > > &existing_nghost_remove,
		vector<vector<vector<int> > > &neigh_cpus)
	{
		// sort vectors for fast search
		vector<vector<size_t> > DuplicatedPoints = tess.GetDuplicatedPoints();
		size_t Nprocs = DuplicatedPoints.size();
		vector<vector<size_t> > sort_indeces(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
		{
			sort_index(DuplicatedPoints[i], sort_indeces[i]);
			sort(DuplicatedPoints[i].begin(), DuplicatedPoints[i].end());
		}
		extensive_remove.clear();
		extensive_remove.resize(Nprocs);
		new_point_remove.clear();
		new_point_remove.resize(Nprocs);
		existing_nghost_remove.clear();
		existing_nghost_remove.resize(Nprocs);
		neigh_cpus.clear();
		neigh_cpus.resize(Nprocs);

		size_t Nremove = ToRemove.size();
		vector<size_t> neigh;
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < Nremove; ++i)
		{
			vector<int> cpus;
			tess.GetNeighbors(ToRemove[i], neigh);
			size_t nneigh = neigh.size();
			// Did we send it to this cpu?
			for (size_t k = 0; k < Nprocs; ++k)
			{
				vector<Vector3D> new_points;
				vector<size_t> exist_neigh;
				vector<size_t>::const_iterator it = binary_find(DuplicatedPoints[k].begin(), DuplicatedPoints[k].end(),
					ToRemove[i]);
				// Is this a boundary point with this cpu?
				if (it != DuplicatedPoints[k].end())
				{
					cpus.push_back(tess.GetDuplicatedProcs()[k]);
					extensive_remove[k].push_back(extensives[ToRemove[i]]);
					for (size_t j = 0; j < nneigh; ++j)
					{
						// Did we send this point already?
						vector<size_t>::const_iterator it2 = binary_find(DuplicatedPoints[k].begin(), DuplicatedPoints[k].end(),
							neigh[j]);
						if (it2 != DuplicatedPoints[k].end())
							exist_neigh.push_back(sort_indeces[k][static_cast<size_t>(it2 - DuplicatedPoints[k].begin())]);
						else
						{
							// New point to send
							if (tess.IsPointOutsideBox(neigh[j]))
								continue;
							if (neigh[j] < Norg)
							{
								new_points.push_back(tess.GetMeshPoint(neigh[j]));
								tess.GetDuplicatedPoints()[k].push_back(neigh[j]);
								DuplicatedPoints[k].push_back(neigh[j]);
								sort(DuplicatedPoints[k].begin(), DuplicatedPoints[k].end());
							}
						}
					}
				}
				new_point_remove[k].push_back(new_points);
				existing_nghost_remove[k].push_back(exist_neigh);
			}
			sort(cpus.begin(), cpus.end());
			cpus = unique(cpus);
			for (size_t k = 0; k < Nprocs; ++k)
				neigh_cpus[k].push_back(cpus);
		}
	}

	vector<size_t> GetMPIRefineSend(Tessellation3D const& tess, vector<size_t> const& newpoints, size_t Ntotal)
	{
		vector<size_t> neigh;
		vector<size_t> to_send;
		size_t Norg = tess.GetPointNo();
		size_t Nrefine = newpoints.size();
		for (size_t i = 0; i < Nrefine; ++i)
		{
			tess.GetNeighbors(newpoints[i], neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] > Norg && neigh[j] < Ntotal && !tess.IsPointOutsideBox(neigh[j]))
				{
					to_send.push_back(newpoints[i]);
					break;
				}
			}
		}
		return to_send;
	}

	vector<vector<size_t> > SendMPIRefine(Tessellation3D const& tess, vector<size_t> const& tosend, vector<vector<Vector3D> > &
		recv_points, vector<vector<vector<int> > > &recv_neigh, vector<vector<size_t> > & splitted_ghosts, vector<size_t> const&
		ToRefine, size_t Ntotal0)
	{
		size_t send_size = tosend.size();
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
		vector<vector<int> > sendDuplicates(Nprocs);
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
						int_temp.push_back(sort_indeces[k][static_cast<int>(it - Nghost[k].begin())]);
					}
				}
				if (!int_temp.empty())
				{
					vector<size_t>::const_iterator it = binary_find(duplicated_points[k].begin(), duplicated_points[k].end(),
						tosend[i]);
					assert(it != duplicated_points[k].end());
					sendDuplicates[k].push_back(static_cast<int>(sort_indeces_duplicated[k]
						[static_cast<size_t>(it - duplicated_points[k].begin())]));
					vector<size_t>::const_iterator it2 = binary_find(ToRefine.begin(), ToRefine.end(), tosend[i]);
					assert(it2 != ToRefine.end());
					sendpoints[k].push_back(tess.GetMeshPoint(static_cast<size_t>(it2 - ToRefine.begin()) + Ntotal0));
					sent_points[k].push_back(tosend[i]);
					sendNghost[k].push_back(int_temp);
				}
			}
			assert(good);
		}
		// Communicate data
		sendDuplicates = MPI_exchange_data(tess.GetDuplicatedProcs(), sendDuplicates);
		splitted_ghosts.clear();
		splitted_ghosts.resize(Nprocs);
		for (size_t i = 0; i < Nprocs; ++i)
			for (size_t j = 0; j < sendDuplicates[i].size(); ++j)
				splitted_ghosts[i].push_back(tess.GetGhostIndeces()[i][static_cast<size_t>(sendDuplicates[i][j])]);
		recv_points = MPI_exchange_data(tess.GetDuplicatedProcs(), sendpoints, tess.GetMeshPoint(0));
		recv_neigh = MPI_exchange_data(tess, sendNghost);
		for (size_t i = 0; i < recv_neigh.size(); ++i)
			for (size_t j = 0; j < recv_neigh[i].size(); ++j)
				for (size_t k = 0; k < recv_neigh[i][j].size(); ++k)
					recv_neigh[i][j][k] = tess.GetDuplicatedPoints()[i][recv_neigh[i][j][k]];
		return sent_points;
	}
#endif

	void RecalcOuterCM(Tessellation3D &tess)
	{
		vector<std::pair<size_t, size_t> > const& allfaceneigh = tess.GetAllFaceNeighbors();
		size_t Nfaces = allfaceneigh.size();
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (allfaceneigh[i].second > Norg && tess.IsPointOutsideBox(allfaceneigh[i].second))
			{
				Vector3D normal = normalize(tess.GetMeshPoint(allfaceneigh[i].second) -
					tess.GetMeshPoint(allfaceneigh[i].first));
				double distance = abs(tess.GetCellCM(allfaceneigh[i].first) - tess.GetFacePoints()
					[tess.GetPointsInFace(i)[0]]);
				tess.GetAllCM()[allfaceneigh[i].second] = tess.GetCellCM(allfaceneigh[i].first) +
					(2 * distance)*normal;
			}
		}
	}

	void CleanOuterPoints(vector<size_t> &neigh, vector<size_t> &nneigh, Tessellation3D const&tess)
	{
		nneigh = RemoveList(nneigh, neigh);
		size_t Norg = tess.GetPointNo();
		vector<size_t> res, res2;
		size_t N = neigh.size();
		res.reserve(N);
		for (size_t i = 0; i < N; ++i)
		{
			if (neigh[i] < Norg || !tess.IsPointOutsideBox(neigh[i]))
				res.push_back(neigh[i]);
		}
		N = nneigh.size();
		res2.reserve(N);
		for (size_t i = 0; i < N; ++i)
		{
			if (nneigh[i] < Norg || !tess.IsPointOutsideBox(nneigh[i]))
				res2.push_back(nneigh[i]);
		}
		neigh = res;
		nneigh = res2;
	}

	void RemoveBadAspectRatio(Tessellation3D const& tess, std::pair<vector<size_t>, vector<Vector3D> > &toremove)
	{
		std::pair<vector<size_t>, vector<Vector3D > > remove_res;
		size_t N = toremove.first.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = tess.GetVolume(toremove.first[i]);
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

	std::pair<vector<size_t>, vector<double> > RemoveBadRemovAspectRatio(vector<double> const& merits, vector<size_t> const& candidates,
		Tessellation3D const& tess)
	{
		vector<size_t> result_names;
		vector<double> result_merits;
		if (merits.size() != candidates.size())
			throw UniversalError("Merits and Candidates don't have same size in RemoveNeighbors");
		vector<size_t> neigh;
		size_t N = merits.size();
		for (size_t i = 0; i < N; ++i)
		{
			bool good = true;
			tess.GetNeighbors(candidates[i], neigh);
			size_t Nn = neigh.size();
			double R = tess.GetWidth(candidates[i]);
			for (size_t j = 0; j < Nn; ++j)
			{
				if (abs(tess.GetMeshPoint(candidates[i]) - tess.GetMeshPoint(neigh[j])) < 0.1*R)
				{
					good = false;
					break;
				}
			}
			if (good)
			{
				result_merits.push_back(merits[i]);
				result_names.push_back(candidates[i]);
			}
		}
		return std::pair<vector<size_t>, vector<double> >(result_names, result_merits);
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
		vector<vector<size_t> > const& duplicated_indeces = tess.GetDuplicatedPoints();
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
							duplicated_indeces[k].end(), neigh[j]);
						if (it != duplicated_indeces[k].end())
						{
							indeces[k].push_back(static_cast<size_t>(it - duplicated_indeces[k].begin()));
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
				sort_index(indeces[i], temp);
				sort(indeces[i].begin(), indeces[i].end());
				all_indeces.insert(all_indeces.end(), indeces[i].begin(), indeces[i].end());
				merit_send[i] = VectorValues(merit_send[i], temp);
				all_merits.insert(all_merits.end(), merit_send[i].begin(), merit_send[i].end());
			}
		}
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
		&nghost_neigh_index, vector<vector<size_t> > &nghost_remove, vector<vector<vector<size_t> > > & duplicate_neigh_index)
	{
		vector<size_t> temp;
		size_t Nprocs = tess.GetDuplicatedProcs().size();
		nghost_neigh_index.resize(Nprocs);
		duplicate_neigh_index.resize(Nprocs);
		vector<vector<Vector3D> > to_add_points(Nprocs), to_add_cm(Nprocs);
		nghost_remove.resize(Nprocs);
		vector<vector<size_t> > duplicated_points = tess.GetDuplicatedPoints();
		vector<vector<size_t> > ghost_points = tess.GetGhostIndeces();
		vector<vector<size_t> > sort_indeces(Nprocs), new_send(Nprocs), sort_indecesg(Nprocs);
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
		}
	}

	void BuildVoronoiMPIRemove(Tessellation3D &tess, size_t &point_to_remove, vector<size_t> &outer_neigh,
		vector<size_t> &duplicate_index, size_t rank_index, Tessellation3D &local, vector<size_t> &nneigh)
	{
		point_to_remove = tess.GetGhostIndeces()[rank_index].at(point_to_remove);
		size_t Nouter = outer_neigh.size();
		for (size_t i = 0; i < Nouter; ++i)
			outer_neigh[i] = tess.GetGhostIndeces()[rank_index].at(outer_neigh[i]);
		Nouter = duplicate_index.size();
		for (size_t i = 0; i < Nouter; ++i)
			duplicate_index[i] = tess.GetDuplicatedPoints()[rank_index].at(duplicate_index[i]);
		sort(duplicate_index.begin(), duplicate_index.end());
		vector<size_t> temp;
		nneigh = outer_neigh;
		for (size_t i = 0; i < duplicate_index.size(); ++i)
		{
			tess.GetNeighbors(duplicate_index[i], temp);
			nneigh.insert(nneigh.end(), temp.begin(), temp.end());
		}
		sort(nneigh.begin(), nneigh.end());
		nneigh = unique(nneigh);
		nneigh = RemoveList(nneigh, duplicate_index);
		RemoveVal(nneigh, point_to_remove);
		vector<Vector3D> points, ghost;
		vector<size_t> toduplicate;
		for (size_t i = 0; i < duplicate_index.size(); ++i)
		{
			points.push_back(tess.GetMeshPoint(duplicate_index[i]));
			toduplicate.push_back(i);
		}
		for (size_t i = 0; i < nneigh.size(); ++i)
			ghost.push_back(tess.GetMeshPoint(nneigh[i]));
		local.BuildNoBox(points, ghost, toduplicate);
	}

	void FixVoronoiRemoveMPI(Tessellation3D &tess, Tessellation3D &local, vector<size_t> &neigh, vector<size_t> &bad_faces,
		vector<double> &dv, vector<size_t> &nneigh, size_t point_to_remove)
	{
		vector<std::pair<size_t, size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		vector<double> & full_area = tess.GetAllArea();
		vector<Vector3D> &full_face_cm = tess.GetAllFaceCM();
		vector<Vector3D>& full_vertices = tess.GetFacePoints();

		size_t Norg = tess.GetPointNo();
		vector<size_t> face_remove;
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			face_remove.clear();
			for (size_t j = 0; j < full_cellfaces[neigh[i]].size(); ++j)
			{
				size_t face = full_cellfaces[neigh[i]][j];
				if (full_faceneigh[face].second == point_to_remove ||
					(full_faceneigh[face].second >= Norg && tess.IsPointOutsideBox(full_faceneigh[face].second)))
				{
					bad_faces.push_back(face);
					face_remove.push_back(face);
				}
			}
			if (!face_remove.empty())
			{
				std::sort(face_remove.begin(), face_remove.end());
				full_cellfaces[neigh[i]] = RemoveList(full_cellfaces[neigh[i]], face_remove);
			}
		}

		// Change CM and volume
		size_t Nneigh = neigh.size();
		dv.clear();
		dv.resize(neigh.size(), 0);
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			dv[i] = local.GetVolume(i) - tess.GetAllVolumes()[neigh[i]];
			tess.GetAllVolumes()[neigh[i]] = local.GetVolume(i);
			tess.GetAllCM()[neigh[i]] = local.GetCellCM(i);
		}
		// Add new faces /Update old
		size_t Nlocal = neigh.size() + nneigh.size() + 4;
		size_t Nfaces = localfaceneigh.size();
		vector<vector<size_t> > neigh_neigh(Nneigh);
		for (size_t i = 0; i < Nneigh; ++i)
			tess.GetNeighbors(neigh[i], neigh_neigh[i]);

		size_t N0, N1;
		size_t nvert = full_vertices.size();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (localfaceneigh[i].first < Nneigh)
			{
				N0 = localfaceneigh[i].first;
				N1 = localfaceneigh[i].second;
				if (N1 < Nneigh || !local.IsPointOutsideBox(N1))
				{
					vector<size_t>::const_iterator it = std::find(neigh_neigh[N0].begin(), neigh_neigh[N0].end(),
						neigh[N1]);
					size_t n0 = neigh[N0];
					size_t n1 = neigh[N1];
					size_t face_index = 0;
					if (it != neigh_neigh[N0].end())
					{
						// We already have this face, just change its points
						face_index = full_cellfaces[n0][static_cast<size_t>(it - neigh_neigh[N0].begin())];
					}
					else
					{
						// New face
						face_index = full_area.size();
						full_area.push_back(0);
						full_facepoints.push_back(vector<size_t>());
						full_face_cm.push_back(Vector3D());
						full_faceneigh.push_back(std::pair<size_t, size_t>(n0, n1));
						full_cellfaces[n0].push_back(face_index);
						if (n1 < Norg)
							full_cellfaces[n1].push_back(face_index);
					}
					full_area[face_index] = local.GetArea(i);
					full_facepoints[face_index] = local.GetPointsInFace(i);
					full_facepoints[face_index] += nvert;
					full_face_cm[face_index] = local.GetAllFaceCM()[i];
					if (n0 > n1)
					{
						size_t ttemp = full_faceneigh[face_index].second;
						full_faceneigh[face_index].second = full_faceneigh[face_index].first;
						full_faceneigh[face_index].first = ttemp;
						FlipVector(full_facepoints[face_index]);
					}
				}
				else
				{
					if (N1 >= Nlocal)
					{
						// We made a new boundary point
						assert(local.IsPointOutsideBox(N1));
						// Add new point
						tess.GetMeshPoints().push_back(local.GetMeshPoint(N1));
						tess.GetAllCM().push_back(local.GetCellCM(N1));
						// Add new face
						full_area.push_back(local.GetArea(i));
						full_facepoints.push_back(local.GetPointsInFace(i));
						full_facepoints.back() += nvert;
						full_face_cm.push_back(local.GetAllFaceCM()[i]);
						size_t n0 = neigh[N0];
						size_t n1 = tess.GetMeshPoints().size() - 1;
						full_cellfaces[n0].push_back(full_area.size() - 1);
						full_faceneigh.push_back(std::pair<size_t, size_t>(n0, n1));
					}
					// else this face should not have changed so no fix is needed
				}
			}
		}
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}
#endif //RICH_MPI

	void CheckCorrect(Tessellation3D const& tess)
	{
		size_t N = tess.GetPointNo();
		for (size_t i = 0; i < N; ++i)
		{
			vector<size_t> const& faces = tess.GetCellFaces(i);
			size_t Nfaces = faces.size();
			double A = 0;
			Vector3D sum;
			for (size_t j = 0; j < Nfaces; ++j)
			{
				double a = tess.GetArea(faces[j]);
				A += a;
				vector<size_t> indeces = tess.GetPointsInFace(faces[j]);
				Vector3D normal2 = tess.GetFaceNeighbors(faces[j]).first == i ? normalize(tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).second) -
					tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).first)) :
					normalize(tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).first) -
						tess.GetMeshPoint(tess.GetFaceNeighbors(faces[j]).second));
				sum += normal2*a;
			}
			vector<size_t> neigh;
			if (abs(sum) > 1e-5*A)
				tess.GetNeighbors(i, neigh);
			assert(abs(sum) < 1e-3*A);
			assert(PointInPoly(tess, tess.GetMeshPoint(i), i));
		}
	}

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

	void BuildLocalVoronoiRemove(Tessellation3D &local, Tessellation3D const& tess, vector<size_t> const& neigh_points,
		vector<size_t> const& nneigh)
	{
		vector<Vector3D> points, ghost;
		vector<size_t> toduplicate;
		for (size_t i = 0; i < neigh_points.size(); ++i)
		{
			points.push_back(tess.GetMeshPoint(neigh_points[i]));
			toduplicate.push_back(i);
		}
		for (size_t i = 0; i < nneigh.size(); ++i)
			ghost.push_back(tess.GetMeshPoint(nneigh[i])); // For MPI this doesn't include Nneigh from other CPU
		local.BuildNoBox(points, ghost, toduplicate);
	}

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

	Conserved3D CalcNewExtensives(Tessellation3D const& tess, Tessellation3D const& local, size_t torefine, vector<size_t> const& neigh,
		vector<ComputationalCell3D> const& cells, EquationOfState const& eos, TracerStickerNames const& tsn,
		vector<Conserved3D> &extensives)
	{
		Conserved3D res;
		PrimitiveToConserved(cells[torefine], tess.GetVolume(torefine), res);
		Conserved3D newpoint(res);
		PrimitiveToConserved(cells[torefine], local.GetVolume(1), extensives[torefine]);
		newpoint -= extensives[torefine];
		assert(newpoint.mass > 0);
		size_t Nreal = neigh.size();
		for (size_t i = 0; i < Nreal; ++i)
		{
			if (tess.GetPointNo() <= neigh[i])
				continue;
			double dV = tess.GetVolume(neigh[i]) - local.GetVolume(i + 2);
			assert(dV > -1e-10*tess.GetVolume(neigh[i]));
			PrimitiveToConserved(cells[neigh[i]], dV, res);
			PrimitiveToConserved(cells[neigh[i]], local.GetVolume(i + 2), extensives[neigh[i]]);
			newpoint += res;
		}
		return newpoint;
	}

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

	void FixVoronoiMPI(Tessellation3D &local, Tessellation3D &tess, size_t refinedghost,
		vector<int> const& refined_neigh, Vector3D const& newpoint, vector<size_t> &bad_faces)
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
					new_face_neigh.second = tess.GetMeshPoints().size();
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
		tess.GetMeshPoints().push_back(newpoint);
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}

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
				if (ghost[i][j] < Ntotal0)
					ghost[i][j] += Nsplit;
		}
	}

	void FixVoronoiRemove(Tessellation3D &local, Tessellation3D &tess, vector<size_t> &neigh, vector<size_t> const& nneigh,
		size_t toremove, vector<double> &dv, vector<size_t> &bad_faces
#ifdef RICH_MPI
		, vector<vector<size_t> > &new_duplicate, vector<vector<size_t> > const& sorted_duplictaed_points
#endif
	)
	{
		vector<std::pair<size_t, size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		vector<double> & full_area = tess.GetAllArea();
		vector<Vector3D> &full_face_cm = tess.GetAllFaceCM();
		vector<Vector3D>& full_vertices = tess.GetFacePoints();
#ifdef RICH_MPI
		vector<size_t> ghost_proc_indeces;
		for (size_t j = 0; j < sorted_duplictaed_points.size(); ++j)
			if (std::binary_search(sorted_duplictaed_points[j].begin(), sorted_duplictaed_points[j].end(), toremove))
				ghost_proc_indeces.push_back(j);
#endif

		size_t Norg = tess.GetPointNo();
		if (toremove < Norg)
		{
			vector<size_t> faces = tess.GetCellFaces(toremove);
			size_t Nfaces = faces.size();
			// Remove old face reference
			for (size_t i = 0; i < Nfaces; ++i)
			{
				size_t other = tess.GetFaceNeighbors(faces[i]).first == toremove ? tess.GetFaceNeighbors(faces[i]).second :
					tess.GetFaceNeighbors(faces[i]).first;
				if (other < Norg)
					RemoveVal(full_cellfaces[other], faces[i]);
				bad_faces.push_back(faces[i]);
			}
			full_cellfaces[toremove].clear();
		}
		vector<size_t> face_remove;
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			face_remove.clear();
			for (size_t j = 0; j < full_cellfaces[neigh[i]].size(); ++j)
			{
				size_t face = full_cellfaces[neigh[i]][j];
				if (full_faceneigh[face].second >= Norg && tess.IsPointOutsideBox(full_faceneigh[face].second))
				{
					bad_faces.push_back(face);
					face_remove.push_back(face);
				}
			}
			if (!face_remove.empty())
			{
				std::sort(face_remove.begin(), face_remove.end());
				full_cellfaces[neigh[i]] = RemoveList(full_cellfaces[neigh[i]], face_remove);
			}
		}

		// Change CM and volume
		size_t Nneigh = neigh.size();
		dv.clear();
		dv.resize(neigh.size(), 0);
		double Vtot = 0;
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			dv[i] = local.GetVolume(i) - tess.GetAllVolumes()[neigh[i]];
			tess.GetAllVolumes()[neigh[i]] = local.GetVolume(i);
			tess.GetAllCM()[neigh[i]] = local.GetCellCM(i);
			Vtot += dv[i];
		}
#ifndef RICH_MPI
		double Vold = tess.GetVolume(toremove);
		assert(Vold > 0.999*Vtot && Vtot > 0.999*Vold);
#endif
		// Add new faces /Update old
		size_t Nlocal = neigh.size() + nneigh.size() + 4;
		size_t Nfaces = localfaceneigh.size();
		vector<vector<size_t> > neigh_neigh(Nneigh);
		for (size_t i = 0; i < Nneigh; ++i)
			tess.GetNeighbors(neigh[i], neigh_neigh[i]);

		size_t N0, N1;
		size_t nvert = full_vertices.size();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (localfaceneigh[i].first < Nneigh)
			{
				N0 = localfaceneigh[i].first;
				N1 = localfaceneigh[i].second;
				if (N1 < Nneigh)
				{
					vector<size_t>::const_iterator it = std::find(neigh_neigh[N0].begin(), neigh_neigh[N0].end(),
						neigh[N1]);
					size_t n0 = neigh[N0];
					size_t n1 = neigh[N1];
					size_t face_index = 0;
					if (it != neigh_neigh[N0].end())
					{
						// We already have this face, just change its points
						face_index = full_cellfaces[n0][static_cast<size_t>(it - neigh_neigh[N0].begin())];
					}
					else
					{
						// New face
						face_index = full_area.size();
						full_area.push_back(0);
						full_facepoints.push_back(vector<size_t>());
						full_face_cm.push_back(Vector3D());
						full_faceneigh.push_back(std::pair<size_t, size_t>(n0, n1));
						full_cellfaces[n0].push_back(face_index);
						if (n1 < Norg)
							full_cellfaces[n1].push_back(face_index);
#ifdef RICH_MPI
						else
							if (!local.IsPointOutsideBox(N1))
							{
								for (size_t k = 0; k < ghost_proc_indeces.size(); ++k)
									new_duplicate[ghost_proc_indeces[k]].push_back(n0);
							}
#endif
					}
					full_area[face_index] = local.GetArea(i);
					full_facepoints[face_index] = local.GetPointsInFace(i);
					full_facepoints[face_index] += nvert;
					full_face_cm[face_index] = local.GetAllFaceCM()[i];
					if (n0 > n1)
					{
						size_t ttemp = full_faceneigh[face_index].second;
						full_faceneigh[face_index].second = full_faceneigh[face_index].first;
						full_faceneigh[face_index].first = ttemp;
						FlipVector(full_facepoints[face_index]);
					}
				}
				else
				{
					if (N1 >= Nlocal)
					{
						// We made a new boundary point
						assert(local.IsPointOutsideBox(N1));
						// Add new point
						tess.GetMeshPoints().push_back(local.GetMeshPoint(N1));
						tess.GetAllCM().push_back(local.GetCellCM(N1));
						// Add new face
						full_area.push_back(local.GetArea(i));
						full_facepoints.push_back(local.GetPointsInFace(i));
						full_facepoints.back() += nvert;
						full_face_cm.push_back(local.GetAllFaceCM()[i]);
						size_t n0 = neigh[N0];
						size_t n1 = tess.GetMeshPoints().size() - 1;
						full_cellfaces[n0].push_back(full_area.size() - 1);
						full_faceneigh.push_back(std::pair<size_t, size_t>(n0, n1));
					}
					// else this face should not have changed so no fix is needed
				}
			}
		}
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}

	void FixExtensiveCellsRemove(vector<ComputationalCell3D> &cells, Tessellation3D const& tess, vector<Conserved3D> &extensives,
		TracerStickerNames const& tsn, vector<double> &dv, size_t toremove, vector<size_t> const& neigh,
		AMRCellUpdater3D const& cu, EquationOfState const& eos
#ifdef RICH_MPI
		, AMRExtensiveUpdater3D const& eu
#endif
	)
	{
		size_t N = dv.size();
		for (size_t i = 0; i < N; ++i)
		{
#ifdef RICH_MPI
			extensives[neigh[i]] += eu.ConvertPrimitveToExtensive3D(cells.at(toremove), eos, dv[i], tsn);
#else
			extensives[neigh[i]] += dv[i] * extensives[toremove] / tess.GetVolume(toremove);
#endif
			cells[neigh[i]] = cu.ConvertExtensiveToPrimitve3D(extensives[neigh[i]], eos, tess.GetVolume(neigh[i]),
				cells[neigh[i]], tsn);
		}
	}
}

AMRCellUpdater3D::~AMRCellUpdater3D(void) {}

AMRExtensiveUpdater3D::~AMRExtensiveUpdater3D(void) {}

Conserved3D SimpleAMRExtensiveUpdater3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
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
	AMRExtensiveUpdater3D* eu) :eos_(eos), refine_(refine), remove_(remove), interp_(slopes), cu_(cu), eu_(eu)
{
	if (!cu)
		cu_ = &scu_;
	if (!eu)
		eu_ = &seu_;
}


void AMR3D::UpdateCellsRefine(Tessellation3D &tess, vector<ComputationalCell3D> &cells, EquationOfState const& eos,
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
	vector<size_t> points_tosend = GetMPIRefineSend(tess, ToRefine.first, Ntotal0);
	vector<vector<size_t> > sent_points = SendMPIRefine(tess, points_tosend, new_points, recv_neigh, splitted_points,
		ToRefine.first, Ntotal0);
	for (size_t i = 0; i < splitted_points.size(); ++i)
	{
		for (size_t j = 0; j < splitted_points[i].size(); ++j)
		{
			tess.GetGhostIndeces()[i].push_back(tess.GetMeshPoints().size());
			BuildLocalVoronoiMPI(vlocal, tess, recv_neigh[i][j], new_points[i][j], splitted_points[i][j]);
			FixVoronoiMPI(vlocal, tess, splitted_points[i][j], recv_neigh[i][j], new_points[i][j], all_bad_faces);
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

void AMR3D::UpdateCellsRemove(Tessellation3D &tess, vector<ComputationalCell3D> &cells, vector<Conserved3D> &extensives,
	EquationOfState const& eos, double time, TracerStickerNames const& tracerstickernames)const
{
	std::pair<vector<size_t>, vector<double> > ToRemove = remove_.ToRemove(tess, cells, time, tracerstickernames);
	ToRemove = RemoveBadRemovAspectRatio(ToRemove.second, ToRemove.first, tess);
	vector<size_t> indeces = sort_index(ToRemove.first);
	ToRemove.second = VectorValues(ToRemove.second, indeces);
	ToRemove.first = VectorValues(ToRemove.first, indeces);
	ToRemove = RemoveNeighbors(ToRemove.second, ToRemove.first, tess);
#ifdef RICH_MPI
	ToRemove = RemoveMPINeighbors(ToRemove.second, ToRemove.first, tess);
	vector<vector<size_t> > sorted_duplicated_points = tess.GetDuplicatedPoints();
	for (size_t i = 0; i < sorted_duplicated_points.size(); ++i)
		sort(sorted_duplicated_points[i].begin(), sorted_duplicated_points[i].end());
	vector<vector<size_t> > new_duplicate(sorted_duplicated_points.size());
#endif
	size_t NRemove = ToRemove.first.size();
	size_t Norg = tess.GetPointNo();
	vector<size_t> neigh, nneigh, temp, bad_faces, temp2;
	vector<double> dv;
	Voronoi3D local(tess.GetBoxCoordinates().first, tess.GetBoxCoordinates().second);
#ifdef debug_amr
	Voronoi3D local2(tess.GetBoxCoordinates().first, tess.GetBoxCoordinates().second);
	vector<Vector3D> mesh = tess.GetMeshPoints();
	mesh.resize(tess.GetPointNo());
	local2.Build(mesh);
#endif

#ifdef RICH_MPI
	vector<vector<vector<size_t> > > nghost_neigh_index, duplicate_neigh_index;
	vector<vector<vector<Vector3D> > > to_add_points;
	vector<vector<size_t> > nghost_remove;
	SendRecvMPIRemoveData(tess, ToRemove.first, nghost_neigh_index, nghost_remove, duplicate_neigh_index);
	for (size_t i = 0; i < nghost_remove.size(); ++i)
	{
		for (size_t j = 0; j < nghost_remove[i].size(); ++j)
		{
			assert(!nghost_neigh_index[i][j].empty() && !duplicate_neigh_index[i][j].empty());
			BuildVoronoiMPIRemove(tess, nghost_remove[i][j], nghost_neigh_index[i][j],
				duplicate_neigh_index[i][j], i, local, nneigh);
			FixVoronoiRemoveMPI(tess, local, duplicate_neigh_index[i][j], bad_faces, dv, nneigh, nghost_remove[i][j]);
			FixExtensiveCellsRemove(cells, tess, extensives, tracerstickernames, dv, nghost_remove[i][j], 
				duplicate_neigh_index[i][j], *cu_, eos,*eu_);
		}
	}

#endif

	for (size_t i = 0; i < NRemove; ++i)
	{
		neigh.clear();
		nneigh.clear();
		tess.GetNeighbors(ToRemove.first[i], temp);
		tess.GetNeighborNeighbors(temp2, ToRemove.first[i]);
		for (size_t j = 0; j < temp.size(); ++j)
		{
			if (temp[j] < Norg)
				neigh.push_back(temp[j]);
			else
				if (!tess.IsPointOutsideBox(temp[j]))
					nneigh.push_back(temp[j]);
		}
		for (size_t j = 0; j < temp2.size(); ++j)
		{
			if (!tess.IsPointOutsideBox(temp2[j]))
				nneigh.push_back(temp2[j]);
		}
		std::sort(neigh.begin(), neigh.end());
		std::sort(nneigh.begin(), nneigh.end());
		BuildLocalVoronoiRemove(local, tess, neigh, nneigh);

#ifdef debug_amr
		vector<double> dv2;
		double vold;
		if (i > 0)
		{
			vold = local2.GetVolume(137 - i);
	}
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			size_t toremove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.begin() + i,
				neigh[j]) - ToRemove.first.begin());
			dv2.push_back(local2.GetVolume(neigh[j] - toremove));
		}
		mesh.erase(mesh.begin() + ToRemove.first[i] - i);
		local2.Build(mesh);
		vector<double> vall, vsmall, vbegin;
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			vbegin.push_back(tess.GetVolume(neigh[j]));
			vsmall.push_back(local.GetVolume(j));
			size_t toremove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.begin() + i + 1,
				neigh[j]) - ToRemove.first.begin());
			vall.push_back(local2.GetVolume(neigh[j] - toremove));
			dv2[j] -= vall.back();
		}
#endif

		FixVoronoiRemove(local, tess, neigh, nneigh, ToRemove.first[i], dv, bad_faces
#ifdef RICH_MPI
			, new_duplicate, sorted_duplicated_points
#endif
		);
		FixExtensiveCellsRemove(cells, tess, extensives, tracerstickernames, dv, ToRemove.first[i], neigh, *cu_, eos
#ifdef RICH_MPI
			, *eu_
#endif
		);
}

	std::sort(bad_faces.begin(), bad_faces.end());
	bad_faces = unique(bad_faces);
	// Remove bad points and faces
	tess.GetPointNo() = Norg - ToRemove.first.size();
	RemoveVector(tess.GetAllArea(), bad_faces);
	RemoveVector(tess.GetAllFaceCM(), bad_faces);
	RemoveVector(tess.GetAllFaceNeighbors(), bad_faces);
	RemoveVector(tess.GetAllPointsInFace(), bad_faces);
	RemoveVector(tess.GetMeshPoints(), ToRemove.first);
	RemoveVector(tess.GetAllCellFaces(), ToRemove.first);
	RemoveVector(tess.GetAllCM(), ToRemove.first);
	RemoveVector(tess.GetAllVolumes(), ToRemove.first);
	RemoveVector(cells, ToRemove.first);
	RemoveVector(extensives, ToRemove.first);
	// Fix face indeces
	vector<std::pair<size_t, size_t> > &face_neigh = tess.GetAllFaceNeighbors();
	size_t Nfaces = face_neigh.size();
	for (size_t i = 0; i < Nfaces; ++i)
	{
		size_t to_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(),
			face_neigh[i].first) - ToRemove.first.begin());
		face_neigh[i].first -= to_remove;
		to_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(),
			face_neigh[i].second) - ToRemove.first.begin());
		face_neigh[i].second -= to_remove;
	}
#ifdef RICH_MPI
	nneigh.clear();
	for (size_t i = 0; i < nghost_remove.size(); ++i)
		for (size_t j = 0; j < nghost_remove[i].size(); ++j)
			nneigh.push_back(nghost_remove[i][j]);
	std::sort(nneigh.begin(), nneigh.end());
	nneigh = unique(nneigh);
	for (size_t i = 0; i < Nfaces; ++i)
	{
		size_t to_remove = static_cast<size_t>(std::lower_bound(nneigh.begin(), nneigh.end(),
			face_neigh[i].first) - nneigh.begin());
		face_neigh[i].first -= to_remove;
		to_remove = static_cast<size_t>(std::lower_bound(nneigh.begin(), nneigh.end(),
			face_neigh[i].second) - nneigh.begin());
		face_neigh[i].second -= to_remove;
	}
#endif
	Norg = tess.GetPointNo();
	vector<vector<size_t> >& allfaces = tess.GetAllCellFaces();
	for (size_t i = 0; i < Norg; ++i)
	{
		size_t Ncell = allfaces[i].size();
		for (size_t j = 0; j < Ncell; ++j)
		{
			size_t to_remove = static_cast<size_t>(std::lower_bound(bad_faces.begin(), bad_faces.end(),
				allfaces[i][j]) - bad_faces.begin());
			allfaces[i][j] -= to_remove;
		}
	}

	// Deal with mpi
#ifdef RICH_MPI
	// Update Nghost
	vector<vector<size_t> > &duplicated_points = tess.GetDuplicatedPoints();
	vector<vector<size_t> > &ghost_indeces = tess.GetGhostIndeces();
	for (size_t i = 0; i < duplicated_points.size(); ++i)
	{
		duplicated_points[i] = RemoveList(duplicated_points[i], ToRemove.first);
		ghost_indeces[i] = RemoveList(ghost_indeces[i], nneigh);
		for (size_t j = 0; j < duplicated_points[i].size(); ++j)
		{
			size_t to_remove = static_cast<size_t>(std::lower_bound(ToRemove.first.begin(), ToRemove.first.end(),
				duplicated_points[i][j]) - ToRemove.first.begin());
			duplicated_points[i][j] -= to_remove;
		}
		for (size_t j = 0; j < ghost_indeces[i].size(); ++j)
		{
			size_t to_remove = static_cast<size_t>(std::lower_bound(nneigh.begin(), nneigh.end(),
				ghost_indeces[i][j]) - nneigh.begin());
			ghost_indeces[i][j] -= to_remove + ToRemove.first.size();
		}
	}
	for (size_t i = 0; i < nneigh.size(); ++i)
		nneigh[i] -= ToRemove.first.size();
	RemoveVector(tess.GetMeshPoints(), nneigh);
	// Update cells and CM
	MPI_exchange_data(tess, tess.GetAllCM(), true);
	MPI_exchange_data(tess, cells, true);
#endif


#ifdef debug_amr
	CheckCorrect(tess);
#endif

}

void AMR3D::operator() (HDSim3D &sim)
{
	UpdateCellsRefine(sim.getTesselation(), sim.getCells(), eos_, sim.getExtensives(), sim.GetTime(),
#ifdef RICH_MPI
		sim.getProcTesselation(),
#endif
		sim.GetTracerStickerNames());
	UpdateCellsRemove(sim.getTesselation(), sim.getCells(), sim.getExtensives(), eos_, sim.GetTime(),
		sim.GetTracerStickerNames());
	// Recalc CM for outerpoints
	RecalcOuterCM(sim.getTesselation());
}

AMR3D::~AMR3D(void) {}