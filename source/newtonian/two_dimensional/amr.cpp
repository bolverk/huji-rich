#include "amr.hpp"
#include "../../tessellation/VoronoiMesh.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif
//#define debug_amr 1

#ifdef debug_amr
#include "hdf5_diagnostics.hpp"
#endif

AMRCellUpdater::~AMRCellUpdater() {}

AMRExtensiveUpdater::~AMRExtensiveUpdater() {}

CellsToRemove::~CellsToRemove() {}

CellsToRefine::~CellsToRefine() {}

AMR::~AMR(void) {}

Extensive SimpleAMRExtensiveUpdater::ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
	double volume, TracerStickerNames const& tracerstickernames) const
{
	Extensive res;
	const double mass = volume*cell.density;
	res.mass = mass;
	size_t N = cell.tracers.size();
	res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell.tracers[i] * mass;
	res.energy = eos.dp2e(cell.density, cell.pressure, cell.tracers,tracerstickernames.tracer_names)*mass +
		0.5*mass*ScalarProd(cell.velocity, cell.velocity);
	res.momentum = mass*cell.velocity;

	return res;
}

SimpleAMRCellUpdater::SimpleAMRCellUpdater(vector<string> toskip) :toskip_(toskip) {}

ComputationalCell SimpleAMRCellUpdater::ConvertExtensiveToPrimitve(const Extensive& extensive, const EquationOfState& eos,
	double volume, ComputationalCell const& old_cell,TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < toskip_.size(); ++i)
		if(safe_retrieve(old_cell.stickers, tracerstickernames.sticker_names,toskip_[i]))
			return old_cell;
	ComputationalCell res;
	const double vol_inv = 1.0 / volume;
	res.density = extensive.mass*vol_inv;
	res.velocity = extensive.momentum / extensive.mass;
	size_t N = extensive.tracers.size();
	res.tracers.resize(N);
	for (size_t i = 0; i < N;++i)
		res.tracers[i]=extensive.tracers[i] / extensive.mass;
	res.stickers = old_cell.stickers;
	res.pressure = eos.de2p(res.density, extensive.energy / extensive.mass - 0.5*ScalarProd(res.velocity, res.velocity),res.tracers,tracerstickernames.tracer_names);
	return res;
}

namespace
{
#ifdef RICH_MPI
	vector<vector<int> > GetSentIndeces(Tessellation const& tess, vector<size_t> const& ToRemove,vector<vector<size_t> > &
		RemoveIndex)
	{
		RemoveIndex.clear();
		vector<vector<int> > sentpoints = tess.GetDuplicatedPoints();
		vector<vector<int> > sort_indeces(sentpoints.size());
		vector<vector<int> > res(sentpoints.size());
		RemoveIndex.resize(sentpoints.size());
		int Nprocs = static_cast<int>(sentpoints.size());
		// sort vectors for fast search
		for (int i = 0; i < Nprocs; ++i)
		{
			sort_index(sentpoints[i], sort_indeces[i]);
			sort(sentpoints[i].begin(), sentpoints[i].end());
		}
		// search the vectors
		int Nremove = static_cast<int>(ToRemove.size());
		for (int i = 0; i < Nremove; ++i)
		{
			for (int j = 0; j < Nprocs; ++j)
			{
				vector<int>::const_iterator it = std::lower_bound(sentpoints[static_cast<size_t>(j)].begin(), sentpoints[static_cast<size_t>(j)].end(),
					ToRemove[static_cast<size_t>(i)]);
				if ((it != sentpoints[static_cast<size_t>(j)].end()) && (*it == static_cast<int>(ToRemove[static_cast<size_t>(i)])))
				{
					res[static_cast<size_t>(j)].push_back(sort_indeces[static_cast<size_t>(j)][static_cast<size_t>(it -
						sentpoints[static_cast<size_t>(j)].begin())]);
					RemoveIndex[static_cast<size_t>(j)].push_back(static_cast<size_t>(i));
				}
			}
		}
		return res;
	}

	void SendRecvOuterMerits(Tessellation const& tess, vector<vector<int> > &sent_indeces,vector<double> const& merits,
		vector<vector<size_t> > const& RemoveIndex,vector<vector<int> > &recv_indeces,vector<vector<double> > &recv_mertis)
	{
		vector<vector<double> > send_merits(sent_indeces.size());
		size_t Nproc = send_merits.size();
		for (size_t i = 0; i < Nproc; ++i)
			send_merits[i]=VectorValues(merits, RemoveIndex[i]);
		recv_indeces = MPI_exchange_data(tess.GetDuplicatedProcs(), sent_indeces);
		recv_mertis = MPI_exchange_data(tess.GetDuplicatedProcs(), send_merits);
	}

	vector<size_t> KeepMPINeighbors(Tessellation const& tess, vector<size_t> const& ToRemove,vector<double> const& merits,
		vector<vector<int> > &recv_indeces,vector<vector<double> > &recv_mertis)
	{
		vector<vector<int> > const& Nghost = tess.GetGhostIndeces();
		size_t Nproc = recv_indeces.size();
		for (size_t i = 0; i < Nproc; ++i)
		{
			recv_indeces[i] = VectorValues(Nghost[i], recv_indeces[i]);
			vector<size_t> temp = sort_index(recv_indeces[i]);
			sort(recv_indeces[i].begin(), recv_indeces[i].end());
			recv_mertis[i] = VectorValues(recv_mertis[i], temp);
		}
		
		int N = tess.GetPointNo();
		size_t Nremove = ToRemove.size();
		vector<int> neigh;
		vector<size_t> RemoveFinal;
		for (size_t i = 0; i < Nremove; ++i)
		{
			bool good = true;
			tess.GetNeighbors(static_cast<int>(ToRemove[i]), neigh);
			size_t Nneigh = neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
			{
				if (neigh[j] >= N)
				{
					for (size_t k = 0; k < Nproc; ++k)
					{
						vector<int>::const_iterator it = binary_find(recv_indeces[k].begin(), recv_indeces[k].end(),neigh[j]);
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

	void ExchangeOuterRemoveData(Tessellation const& tess, vector<size_t> const& ToRemove,vector<vector<int> >
		&to_check,vector<vector<Vector2D> > &chulls_res,vector<Extensive> const& extensives,vector<Extensive> & mpi_extensives,
		CacheData const& cd)
	{
		chulls_res.clear();
		to_check.clear();
		mpi_extensives.clear();
		vector<vector<int> > Nghost = tess.GetGhostIndeces();
		vector<vector<int> > sort_indeces(Nghost.size());
		int Nprocs = static_cast<int>(Nghost.size());
		// sort vectors for fast search
		for (int i = 0; i < Nprocs; ++i)
		{
			sort_index(Nghost[i], sort_indeces[i]);
			sort(Nghost[i].begin(), Nghost[i].end());
		}
		vector<vector<Extensive> > myremove(Nprocs);
		vector<vector<vector<Vector2D> > > chulls(Nprocs);
		vector<vector<vector<int> > > remove_neigh(Nprocs);
		size_t Nremove = ToRemove.size();
		int Np = tess.GetPointNo();
		vector<int> neigh;
		vector<Vector2D> tempV2D;
		for (size_t i = 0; i < Nremove; ++i)
		{
			bool calculated = false;
			tess.GetNeighbors(static_cast<int>(ToRemove[i]), neigh);
			size_t Nneigh = neigh.size();
			for (size_t k = 0; k < static_cast<size_t>(Nprocs); ++k)
			{
				vector<int> temp_add;
				bool should_add=false;
				for (size_t j = 0; j < Nneigh; ++j)
				{
					if (neigh[j] >= Np && tess.GetOriginalIndex(neigh[j]) != static_cast<int>(ToRemove[i]))
					{
						vector<int>::const_iterator it = binary_find(Nghost[k].begin(), Nghost[k].end(), neigh[j]);
						if (it != Nghost[k].end())
						{
							should_add=true;
							if (!calculated)
							{
								ConvexHull(tempV2D, tess, static_cast<int>(ToRemove[i]));
								calculated = true;
							}
							temp_add.push_back(sort_indeces[k][static_cast<size_t>(it - Nghost[k].begin())]);
						}
					}
				}
				if (should_add)
				{
					chulls[k].push_back(tempV2D);
					remove_neigh[k].push_back(temp_add);
					myremove[k].push_back((1.0/cd.volumes[ToRemove[i]]) * extensives[ToRemove[i]]);
				}
			}
		}
		// Exchange the data
		chulls=MPI_exchange_data(tess, chulls, tess.GetMeshPoint(0));
		remove_neigh=MPI_exchange_data(tess, remove_neigh);
		vector<vector<Extensive> > mpi_recv_extensives = MPI_exchange_data(tess.GetDuplicatedProcs(), myremove,extensives[0]);
		mpi_extensives = CombineVectors(mpi_recv_extensives);
		//convert  remove_neigh to indeces of real points via duplciated points
		vector<vector<int> > duplicated_points = tess.GetDuplicatedPoints();
		// sort vectors for fast search
		/*for (int i = 0; i < Nprocs; ++i)
		{
			sort_index(duplicated_points[i], sort_indeces[i]);
			sort(duplicated_points[i].begin(), duplicated_points[i].end());
		}*/
		for (size_t i = 0; i < static_cast<size_t>(Nprocs); ++i)
		{
			for (size_t j = 0; j < remove_neigh[i].size(); ++j)
			{
				for (size_t k = 0; k < remove_neigh[i][j].size(); ++k)
				{
					/*vector<int>::const_iterator it = binary_find(duplicated_points[i].begin(), duplicated_points[i].end(), remove_neigh[i][j][k]);
					if (it != duplicated_points[i].end())
					{
					size_t loc = static_cast<size_t>(it - duplicated_points[i].begin());
						remove_neigh[i][j][k] = duplicated_points[i][sort_indeces[i][loc]];
					}*/
					remove_neigh[i][j][k]=duplicated_points[i][remove_neigh[i][j][k]];
				}
			}
		}
		to_check = CombineVectors(remove_neigh);
		chulls_res = CombineVectors(chulls);
	}

	void DealWithMPINeighbors(Tessellation const& tess, vector<size_t> &ToRemove, vector<double> &merits,
		vector<vector<int> > &to_check, vector<vector<Vector2D> > &chulls_res,
		vector<Extensive> const& extensives, vector<Extensive> & mpi_extensives,CacheData const& cd)
	{
		vector<vector<size_t> > RemoveIndex;
		vector<vector<int> > sent_indeces = GetSentIndeces(tess, ToRemove, RemoveIndex);
		vector<vector<int> > recv_indeces;
		vector<vector<double> > recv_mertis;
		SendRecvOuterMerits(tess, sent_indeces, merits, RemoveIndex, recv_indeces, recv_mertis);
		ToRemove = KeepMPINeighbors(tess, ToRemove, merits, recv_indeces, recv_mertis);
		ExchangeOuterRemoveData(tess, ToRemove, to_check, chulls_res,extensives, mpi_extensives,cd);
	}
#endif

	double AreaOverlap(vector<Vector2D> const& poly0, vector<Vector2D> const& poly1,PhysicalGeometry const& pg)
	{
		// returns the overlap of p1 with p0
		using namespace ClipperLib;
		Paths subj(1), clip(1), solution;
		double maxi = 0;
		Vector2D cm0;
		for (size_t i = 0; i < poly0.size(); ++i)
			cm0 += poly0[i];
		cm0 = cm0/static_cast<double>(poly0.size());
		for (size_t i = 0; i < poly0.size(); ++i)
			maxi = std::max(maxi, std::max(std::abs(poly0[i].x-cm0.x), std::abs(poly0[i].y-cm0.y)));
		for (size_t i = 0; i < poly1.size(); ++i)
			maxi = std::max(maxi, std::max(std::abs(poly1[i].x-cm0.x), std::abs(poly1[i].y-cm0.y)));
		int maxscale = static_cast<int>(log10(maxi) + 10);

		subj[0].resize(poly0.size());
		clip[0].resize(poly1.size());
		for (size_t i = 0; i < poly0.size(); ++i)
			subj[0][i] = IntPoint(static_cast<cInt>((poly0[i].x-cm0.x)*pow(10.0, 18 - maxscale)), static_cast<cInt>((poly0[i].y-cm0.y)*pow(10.0, 18 - maxscale)));
		for (size_t i = 0; i < poly1.size(); ++i)
			clip[0][i] = IntPoint(static_cast<cInt>((poly1[i].x-cm0.x)*pow(10.0, 18 - maxscale)), static_cast<cInt>((poly1[i].y-cm0.y)*pow(10.0, 18 - maxscale)));

		//perform intersection ...
		Clipper c;
		c.AddPaths(subj, ptSubject, true);
		c.AddPaths(clip, ptClip, true);
		c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
		if (!solution.empty())
		{
			if (solution[0].size() > 2)
			{
				vector<Vector2D> inter(solution[0].size());
				double factor = pow(10.0, maxscale - 18);
				for (size_t i = 0; i < solution[0].size(); ++i)
					inter[i].Set(static_cast<double>(solution[0][i].X) * factor + cm0.x, static_cast<double>(solution[0][i].Y) * factor + cm0.y);
				return pg.calcVolume(inter);
			}
		}
		return 0;
	}

	Vector2D FixInDomain(OuterBoundary const& obc, Vector2D &point)
	{
		Vector2D res;
		if (point.x > obc.GetGridBoundary(Right))
		{
			point.x -= obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
			res.x += obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
		}
		if (point.x < obc.GetGridBoundary(Left))
		{
			point.x += obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
			res.x -= obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
		}
		if (point.y > obc.GetGridBoundary(Up))
		{
			point.y -= obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
			res.y += obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
		}
		if (point.y < obc.GetGridBoundary(Down))
		{
			point.y += obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
			res.y -= obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
		}
		return res;
	}
}

namespace
{
	double GetAspectRatio(Tessellation const& tess, int index)
	{
		vector<int> const& edges = tess.GetCellEdges(index);
		double L = 0;
		for (size_t i = 0; i < edges.size(); ++i)
			L += tess.GetEdge(edges[i]).GetLength();
		return 4 * 3.14*tess.GetVolume(index) / (L*L);
	}

	void RemoveNeighbors(vector<double> &merits, vector<size_t> &candidates, Tessellation const& tess)
	{
		vector<size_t> candidates_new;
		vector<double> merits_new;
		if (merits.size() != candidates.size())
			throw UniversalError("Merits and Candidates don't have same size in RemoveNeighbors");
		// Make sure there are no neighbors
		vector<size_t> bad_neigh;
		vector<int> neigh;
		for (size_t i = 0; i < merits.size(); ++i)
		{
			bool good = true;
			tess.GetNeighbors(static_cast<int>(candidates[i]),neigh);
			size_t nneigh = neigh.size();
			if (find(bad_neigh.begin(), bad_neigh.end(), candidates[i]) !=
				bad_neigh.end())
				good = false;
			else
			{
				for (size_t j = 0; j < nneigh; ++j)
				{
					if (binary_search(candidates.begin(), candidates.end(), tess.GetOriginalIndex(neigh[j])))
					{
						if (merits[i] < merits[static_cast<size_t>(lower_bound(candidates.begin(), candidates.end(),
							tess.GetOriginalIndex(neigh[j])) - candidates.begin())])
						{
							good = false;
							break;
						}
						if (fabs(merits[i] - merits[static_cast<size_t>(lower_bound(candidates.begin(), candidates.end(),
							tess.GetOriginalIndex(neigh[j])) - candidates.begin())]) < 1e-9)
						{
							if (find(bad_neigh.begin(), bad_neigh.end(), tess.GetOriginalIndex(neigh[j])) == bad_neigh.end())
								bad_neigh.push_back(static_cast<size_t>(tess.GetOriginalIndex(neigh[j])));
						}
					}
				}
			}
			if (good)
			{
				candidates_new.push_back(candidates[i]);
				merits_new.push_back(merits[i]);
			}
		}
		candidates = candidates_new;
		merits = merits_new;
	}

	Extensive GetNewExtensive(vector<Extensive> const& extensives,Tessellation const& tess,size_t N,size_t location,
		vector<Vector2D> const& moved,vector<int> const& real_neigh,vector<vector<Vector2D> >  const& Chull,
		vector<ComputationalCell> const& cells,EquationOfState const& eos, AMRExtensiveUpdater const& eu,
		double &TotalVolume,Tessellation const& oldtess,bool periodic,TracerStickerNames const& tracerstickernames,PhysicalGeometry const& pg,
		CacheData const& cd)
	{
		Extensive NewExtensive(extensives[0].tracers);
		vector<Vector2D> temp;
		ConvexHull(temp, tess, static_cast<int>(N + location));
		temp = temp + moved[location];
		const double vv = cd.volumes[N + location];
		stack<int> tocheck;
		stack<vector<Vector2D> > tocheck_hull;
		vector<int> checked(real_neigh);
		vector<int> neightemp;
		for (size_t j = 0; j < Chull.size(); ++j)
		{
			double v = AreaOverlap(temp, Chull[j],pg);
			if (v > vv*1e-8)
			{
				NewExtensive += eu.ConvertPrimitveToExtensive(cells[static_cast<size_t>(real_neigh[j])], eos, v,tracerstickernames);
				TotalVolume += v;
#ifdef RICH_MPI
				if (oldtess.GetOriginalIndex(real_neigh[j]) > tess.GetPointNo() || oldtess.GetOriginalIndex(real_neigh[j])==real_neigh[j])
					continue;
#endif
				oldtess.GetNeighbors(real_neigh[j],neightemp);
				for (size_t k = 0; k < neightemp.size(); ++k)
				{
#ifdef RICH_MPI
					if (oldtess.GetOriginalIndex(neightemp[k]) > tess.GetPointNo() || oldtess.GetOriginalIndex(neightemp[k]) == neightemp[k])
						continue;
#endif
					if (std::find(checked.begin(), checked.end(), oldtess.GetOriginalIndex(neightemp[k])) == checked.end())
					{
						vector<Vector2D> temp2;
						ConvexHull(temp2, oldtess, oldtess.GetOriginalIndex(neightemp[k]));
						tocheck.push(oldtess.GetOriginalIndex(neightemp[k]));
						tocheck_hull.push(temp2);
						checked.push_back(oldtess.GetOriginalIndex(neightemp[k]));
					}
				}
			}		
		}
		while (!tocheck.empty())
		{
			vector<Vector2D> chull = tocheck_hull.top();
			int cur_check=tocheck.top();
			tocheck.pop();
			tocheck_hull.pop();
			double v = AreaOverlap(temp, chull,pg);
			if (v > vv*1e-8)
			{
				NewExtensive += eu.ConvertPrimitveToExtensive(cells[static_cast<size_t>(cur_check)], eos, v,tracerstickernames);
				TotalVolume += v;
				oldtess.GetNeighbors(cur_check,neightemp);
				for (size_t k = 0; k < neightemp.size(); ++k)
				{
#ifdef RICH_MPI
					if (oldtess.GetOriginalIndex(neightemp[k]) > tess.GetPointNo())
						continue;
#endif
					if (std::find(checked.begin(), checked.end(), oldtess.GetOriginalIndex(neightemp[k])) == checked.end())
					{
						vector<Vector2D> temp2;
						ConvexHull(temp2, oldtess, oldtess.GetOriginalIndex(neightemp[k]));
						tocheck.push(oldtess.GetOriginalIndex(neightemp[k]));
						tocheck_hull.push(temp2);
						checked.push_back(oldtess.GetOriginalIndex(neightemp[k]));
					}
				}
			}
		}
		const double eps = periodic ? 1e-2 : 1e-5;
		if (vv > (1 + eps)*TotalVolume || vv < (1 - eps)*TotalVolume)
		{
			std::cout << "In refine Real volume: " << vv << " AMR volume: " << TotalVolume << std::endl;
#ifndef RICH_MPI
			UniversalError eo("Not same volume in amr refine");
			eo.AddEntry("location",static_cast<double>(location));
			throw eo;
#endif
		}
		return NewExtensive;
	}

	void GetToCheck(Tessellation const& tess, int ToRefine, vector<int> &real_neigh, vector<vector<Vector2D> > &Chull)
	{
		vector<int> neigh = tess.GetNeighbors(ToRefine);
		vector<Vector2D> temp;
		ConvexHull(temp, tess, static_cast<int>(ToRefine));
		Chull.push_back(temp);
		real_neigh.push_back(static_cast<int>(ToRefine));
		size_t N = static_cast<size_t>(tess.GetPointNo());
		vector<Vector2D> moved;
		moved.push_back(Vector2D(0, 0));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			if (neigh[j] < static_cast<int>(N))
			{
				real_neigh.push_back(neigh[j]);
				ConvexHull(temp, tess, neigh[j]);
				Chull.push_back(temp);
				moved.push_back(Vector2D(0, 0));
			}
			else
			{
#ifndef RICH_MPI
				// Is it not a rigid wall?
				if (tess.GetOriginalIndex(neigh[j]) != static_cast<int>(ToRefine))
				{
					real_neigh.push_back(tess.GetOriginalIndex(neigh[j]));
					ConvexHull(temp, tess, real_neigh.back());
					temp = temp + (tess.GetMeshPoint(neigh[j]) -
						tess.GetMeshPoint(tess.GetOriginalIndex(neigh[j])));
					Chull.push_back(temp);
					moved.push_back((tess.GetMeshPoint(neigh[j]) - tess.GetMeshPoint(tess.GetOriginalIndex(neigh[j]))));
				}
#endif
			}
		}
		// Get neighbor neighbors
		const size_t NN = real_neigh.size();
		vector<int> neigh_temp;
		for (size_t j = 0; j < NN; ++j)
		{
			tess.GetNeighbors(real_neigh[j],neigh_temp);
			for (size_t i = 0; i < neigh_temp.size(); ++i)
			{
				if (std::find(real_neigh.begin(), real_neigh.end(), tess.GetOriginalIndex(neigh_temp[i])) == real_neigh.end())
				{
					if (neigh_temp[i] < static_cast<int>(N))
					{
						real_neigh.push_back(neigh_temp[i]);
						ConvexHull(temp, tess, neigh_temp[i]);
						Chull.push_back(temp+moved[j]);
					}
					else
					{
#ifndef RICH_MPI
						// Is it not a rigid wall?
						if (tess.GetOriginalIndex(neigh_temp[i]) != static_cast<int>(real_neigh[j]))
						{
							real_neigh.push_back(tess.GetOriginalIndex(neigh_temp[i]));
							ConvexHull(temp, tess, real_neigh.back());
							temp = temp + (tess.GetMeshPoint(neigh_temp[i]) -
								tess.GetMeshPoint(tess.GetOriginalIndex(neigh_temp[i]))) + moved[j];
							Chull.push_back(temp);
						}
#endif
					}
				}
			}
		}
	}

	void AddCurrentNeigh(Tessellation const& tess,
		size_t location, vector<vector<Vector2D> > &Chull, vector<int> &real_neigh,vector<Vector2D> const& Moved,
		Tessellation const& oldtess)
	{
		int oldpointnumber = oldtess.GetPointNo();
		int newindex = static_cast<int>(location) + oldpointnumber;
		vector<int> neigh = tess.GetNeighbors(newindex);
		vector<int> neighcopy(real_neigh);
		sort(neighcopy.begin(), neighcopy.end());
		vector<Vector2D> temp;
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			// Is it a new point? Should deal with mpi as well
			if (tess.GetOriginalIndex(neigh[i]) >= oldpointnumber)
				continue;
			// Rigid wall?
			if (tess.GetOriginalIndex(neigh[i]) == newindex)
				continue;
			// Do we already have this point?
			if (binary_search(neighcopy.begin(), neighcopy.end(), tess.GetOriginalIndex(neigh[i])))
				continue;
			real_neigh.push_back(tess.GetOriginalIndex(neigh[i]));
			ConvexHull(temp, oldtess, real_neigh.back());
			// Periodic?
			if (neigh[i] > tess.GetPointNo())
				temp = temp + (tess.GetMeshPoint(neigh[i]) -
					tess.GetMeshPoint(tess.GetOriginalIndex(neigh[i])))+Moved[location];
			Chull.push_back(temp);
		}
	}

	void ConservedSingleCell(Tessellation const& oldtess, Tessellation const& tess, size_t ToRefine, size_t &location,
		LinearGaussImproved *interp, EquationOfState const& eos, vector<ComputationalCell> &cells,
		vector<pair<size_t, Vector2D> > const& NewPoints, vector<Extensive> const& extensives,
		AMRExtensiveUpdater const& eu, AMRCellUpdater const& cu,vector<Vector2D> const& moved,bool periodic,
		TracerStickerNames const& tracerstickernames,PhysicalGeometry const& pg,CacheData const& cd)
	{
		vector<int> real_neigh;
		vector<vector<Vector2D> > Chull;
		GetToCheck(oldtess, static_cast<int>(ToRefine),real_neigh,Chull);
		size_t N = static_cast<size_t>(oldtess.GetPointNo());
		while (location < NewPoints.size() && NewPoints[location].first == ToRefine)
		{
			AddCurrentNeigh(tess,location, Chull, real_neigh, moved, oldtess);
			double TotalVolume=0;
			Extensive NewExtensive = GetNewExtensive(extensives, tess, N, location, moved, real_neigh, Chull,
				cells, eos, eu,TotalVolume,oldtess,periodic,tracerstickernames,pg,cd);
			cells.push_back(cu.ConvertExtensiveToPrimitve(NewExtensive, eos, TotalVolume, cells[ToRefine],tracerstickernames));
			if (interp != 0)
				interp->GetSlopesUnlimited().push_back(interp->GetSlopesUnlimited()[ToRefine]);
			++location;
		}

	}

}


void AMR::GetNewPoints2(vector<size_t> const& ToRefine, Tessellation const& tess,
	vector<std::pair<size_t, Vector2D> > &NewPoints, vector<Vector2D> &Moved,
	OuterBoundary const& obc)const
{
	const double small = 1e-3;
	NewPoints.clear();
	NewPoints.reserve(ToRefine.size());
	Moved.clear();
	Moved.reserve(ToRefine.size());
	vector<int> neigh;
	for (size_t i = 0; i < ToRefine.size(); ++i)
	{
		// Split only if cell is rather round
		const double R = tess.GetWidth(static_cast<int>(ToRefine[i]));
		Vector2D const& myCM = tess.GetCellCM(static_cast<int>(ToRefine[i]));
		Vector2D const& mypoint = tess.GetMeshPoint(static_cast<int>(ToRefine[i]));
		if (GetAspectRatio(tess, static_cast<int>(ToRefine[i])) < 0.65)
			continue;
		if (myCM.distance(mypoint) > 0.2*R)
			continue;
		tess.GetNeighbors(static_cast<int>(ToRefine[i]),neigh);
		double otherdist = mypoint.distance(tess.GetMeshPoint(neigh[0]));
		size_t max_ind = 0;
		for (size_t j = 1; j < neigh.size(); ++j)
		{
			Vector2D const& otherpoint = tess.GetMeshPoint(neigh[j]);
			if (otherpoint.distance(mypoint) > otherdist)
			{
				max_ind = j;
				otherdist = otherpoint.distance(mypoint);
			}
		}
		Vector2D direction = tess.GetMeshPoint(neigh[max_ind]) - mypoint;
		Vector2D candidate = mypoint+small*direction/abs(direction);
		Moved.push_back(FixInDomain(obc, candidate));
		NewPoints.push_back(std::pair<size_t, Vector2D>(ToRefine[i], candidate));
	}
}


void AMR::GetNewPoints(vector<size_t> const& ToRefine, Tessellation const& tess,
	vector<std::pair<size_t, Vector2D> > &NewPoints, vector<Vector2D> &Moved,
	OuterBoundary const& obc
#ifdef RICH_MPI
	, vector<Vector2D> const& proc_chull
#endif
	)const
{
	// ToRefine should be sorted
	size_t N = static_cast<size_t>(tess.GetPointNo());
	NewPoints.clear();
	NewPoints.reserve(ToRefine.size() * 7);
	Moved.clear();
	Moved.reserve(ToRefine.size() * 7);
	vector<int> neigh;
	for (size_t i = 0; i < ToRefine.size(); ++i)
	{
		// Split only if cell is rather round
		const double R = tess.GetWidth(static_cast<int>(ToRefine[i]));
		Vector2D const& mypoint = tess.GetCellCM(static_cast<int>(ToRefine[i]));
		if (GetAspectRatio(tess, static_cast<int>(ToRefine[i])) < 0.65)
			continue;

		tess.GetNeighbors(static_cast<int>(ToRefine[i]),neigh);
		vector<int> const& edges = tess.GetCellEdges(static_cast<int>(ToRefine[i]));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			Vector2D const& otherpoint = tess.GetCellCM(neigh[j]);
			if (otherpoint.distance(mypoint) < 1.75*R)
				continue;
			if (tess.GetEdge(edges[j]).GetLength() < 0.5*R)
				continue;
			Vector2D candidate = 0.75*mypoint + 0.25*otherpoint;
			if (candidate.distance(tess.GetMeshPoint(static_cast<int>(ToRefine[i]))) < 0.5*R)
				continue;
			if (neigh[j] < static_cast<int>(N) && candidate.distance(tess.GetMeshPoint(neigh[j])) < tess.GetWidth(neigh[j])*0.5)
				continue;
			// Make sure not to split neighboring cells, this causes bad aspect ratio
			if (std::binary_search(ToRefine.begin(), ToRefine.end(), static_cast<size_t>(neigh[j])))
			{
				if (static_cast<size_t>(neigh[j]) > ToRefine[i])
				{
					candidate = 0.5*mypoint + 0.5*otherpoint;
					if (candidate.distance(tess.GetMeshPoint(static_cast<int>(ToRefine[i]))) < 0.5*R ||
						(neigh[j] < static_cast<int>(N) && candidate.distance(tess.GetMeshPoint(neigh[j])) < tess.GetWidth(neigh[j])*0.5))
					{
						candidate = 0.5*tess.GetMeshPoint(static_cast<int>(ToRefine[i])) + 0.5*tess.GetMeshPoint(neigh[j]);
						if (candidate.distance(tess.GetMeshPoint(static_cast<int>(ToRefine[i]))) < 0.5*R ||
							(neigh[j] < static_cast<int>(N) && candidate.distance(tess.GetMeshPoint(neigh[j])) < tess.GetWidth(neigh[j])*0.5))
							continue;
					}
				}
				else
					continue;
			}
#ifdef RICH_MPI
			if (!PointInCell(proc_chull, candidate))
				continue;
#endif
			Moved.push_back(FixInDomain(obc, candidate));
			NewPoints.push_back(std::pair<size_t, Vector2D>(ToRefine[i], candidate));
		}
	}
}

ConservativeAMR::ConservativeAMR
(CellsToRefine const& refine,
	CellsToRemove const& remove,
	bool periodic,
	LinearGaussImproved *slopes,
	AMRCellUpdater* cu,
	AMRExtensiveUpdater* eu) :
	refine_(refine),
	remove_(remove),
	scu_(vector<string>()),
	seu_(),
	periodic_(periodic),
	cu_(cu),
	eu_(eu),
	interp_(slopes)
{
	if (!cu) 
		cu_ = &scu_;
	if(!eu)
		eu_ = &seu_;
}

void ConservativeAMR::UpdateCellsRefine
(Tessellation &tess,
	OuterBoundary const& obc,
	vector<ComputationalCell> &cells,
	EquationOfState const& eos,
	vector<Extensive> &extensives,
	double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames,
	CacheData const& cd,
	PhysicalGeometry const& pg)const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	// Find the primitive of each point and the new location
	vector<std::pair<size_t, Vector2D> > NewPoints;
	vector<Vector2D> Moved;
	vector<size_t> ToRefine = refine_.ToRefine(tess, cells, time,tracerstickernames);
#ifndef RICH_MPI
	if (ToRefine.empty())
		return;
#endif // RICH_MPI
	sort(ToRefine.begin(), ToRefine.end());
	ToRefine=unique(ToRefine);
#ifdef RICH_MPI
	vector<double> merit;
	ToRefine = RemoveNearBoundaryPoints(ToRefine, tess,merit);
	vector<Vector2D> chull;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ConvexHull(chull, proctess, rank);
#endif
	GetNewPoints(ToRefine, tess, NewPoints, Moved, obc
#ifdef RICH_MPI
		, chull
#endif
		);
	// save copy of old tessellation
	boost::scoped_ptr<Tessellation> oldtess(tess.clone());

	// Rebuild tessellation
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);

	for (size_t i = 0; i < NewPoints.size(); ++i)
		cor.push_back(NewPoints[i].second);
#ifdef debug_amr
	WriteTess(tess, "c:/vold.h5");
#endif

#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif
	cd.reset();
#ifdef debug_amr
	WriteTess(tess, "c:/vnew.h5");
#endif


	size_t location = 0;
	for (size_t i = 0; i < ToRefine.size(); ++i)
		ConservedSingleCell(*oldtess, tess, ToRefine[i], location, interp_, eos, cells, NewPoints, extensives, *eu_,
			*cu_,Moved,periodic_,tracerstickernames,pg,cd);
	extensives.resize(N + NewPoints.size());
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, cd.volumes[i],tracerstickernames);
#ifdef RICH_MPI
	MPI_exchange_data(tess, cells, true);
#endif
}

void ConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames,
	CacheData const& cd,
	PhysicalGeometry const& pg)const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	// Rebuild tessellation
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);
	std::pair<vector<size_t>, vector<double> > ToRemovepair = remove_.ToRemove(tess, cells, time,tracerstickernames);
#ifndef RICH_MPI
	if (ToRemovepair.first.empty())
		return;
#endif // RICH_MPI

	// Clean up vectors
	vector<int> indeces;
	sort_index(ToRemovepair.first,indeces);
	sort(ToRemovepair.first.begin(), ToRemovepair.first.end());
	ToRemovepair.second = VectorValues(ToRemovepair.second, indeces);
	indeces = unique_index(ToRemovepair.first);
	ToRemovepair.first = unique(ToRemovepair.first);
	ToRemovepair.second = VectorValues(ToRemovepair.second, indeces);
	RemoveNeighbors(ToRemovepair.second, ToRemovepair.first, tess);
#ifdef RICH_MPI
	vector<vector<int> > mpi_check;
	vector<vector<Vector2D> > chulls_mpi;
	vector<Extensive> mpi_extensives;
	DealWithMPINeighbors(tess, ToRemovepair.first, ToRemovepair.second, mpi_check, chulls_mpi,extensives, mpi_extensives,cd);
#endif

	// save copy of old tessellation
	boost::scoped_ptr<Tessellation> oldtess(tess.clone());
	vector<double> old_vol(ToRemovepair.first.size());
	for (size_t i = 0; i < ToRemovepair.first.size(); ++i)
		old_vol[i] = cd.volumes[ToRemovepair.first[i]];
	RemoveVector(cor, ToRemovepair.first);

#ifdef debug_amr
	WriteTess(tess, "c:/vold.h5");
#endif

#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif
	cd.reset();
#ifdef debug_amr
	WriteTess(tess, "c:/vnew.h5");
#endif

	// Fix the extensives
	vector<Vector2D> temp;
	vector<Vector2D> chull;
	for (size_t i = 0; i < ToRemovepair.first.size(); ++i)
	{
		vector<int> neigh = oldtess->GetNeighbors(static_cast<int>(ToRemovepair.first[i]));
		ConvexHull(chull, *oldtess, static_cast<int>(ToRemovepair.first[i]));
		const double TotalV = old_vol[i];
		temp.clear();
		double dv = 0;
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			size_t toadd = static_cast<size_t>(lower_bound(ToRemovepair.first.begin(), ToRemovepair.first.end(),
				oldtess->GetOriginalIndex(neigh[j])) - ToRemovepair.first.begin());
			temp.clear();
			if (neigh[j] < static_cast<int>(N))
			{
				ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(neigh[j]) - toadd));
			}
			else
			{
				if (oldtess->GetOriginalIndex(neigh[j]) != static_cast<int>(ToRemovepair.first[i]))
				{
#ifdef RICH_MPI
					continue;
#endif
					ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(
						oldtess->GetOriginalIndex(neigh[j]) )- toadd));
					temp = temp + (oldtess->GetMeshPoint(neigh[j]) -
						oldtess->GetMeshPoint(oldtess->GetOriginalIndex(neigh[j])));
				}
			}
			if (!temp.empty())
			{
				const double v = AreaOverlap(chull, temp,pg);
				dv += v;
				extensives[static_cast<size_t>(oldtess->GetOriginalIndex(neigh[j]))] += (v / TotalV)*extensives[ToRemovepair.first[i]];
			}
		}
		/*if (dv > (1 + 1e-5)*TotalV || dv < (1 - 1e-5)*TotalV)
		{
			std::cout << "Real volume: " << TotalV << " AMR volume: " << dv << std::endl;
			throw UniversalError("Not same volume in amr remove");
		}*/
	}
#ifdef RICH_MPI
	for (size_t i = 0; i < mpi_check.size(); ++i)
	{
		for (size_t j = 0; j < mpi_check[i].size(); ++j)
		{
			size_t toadd = static_cast<size_t>(lower_bound(ToRemovepair.first.begin(), ToRemovepair.first.end(),
				oldtess->GetOriginalIndex(mpi_check[i][j])) - ToRemovepair.first.begin());
			ConvexHull(chull, tess, mpi_check[i][j]-static_cast<int>(toadd));
			const double v = AreaOverlap(chull, chulls_mpi[i], pg);
			extensives[static_cast<size_t>(mpi_check[i][j])] += v*mpi_extensives[i];
		}
	}
#endif

	RemoveVector(extensives, ToRemovepair.first);
	RemoveVector(cells, ToRemovepair.first);

	N = static_cast<size_t>(tess.GetPointNo());
	cells.resize(static_cast<size_t>(tess.GetPointNo()));
	for (size_t i = 0; i < N; ++i)
		cells[i] = cu_->ConvertExtensiveToPrimitve(extensives[i], eos,cd.volumes[i], cells[i],tracerstickernames);
#ifdef RICH_MPI
	MPI_exchange_data(tess, cells, true);
#endif
}

void ConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames(),sim.getCacheData(),sim.getPhysicalGeometry());
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames(), sim.getCacheData(), sim.getPhysicalGeometry());
	// redo cache data
	sim.getCacheData().reset();
}



NonConservativeAMR::NonConservativeAMR
(CellsToRefine const& refine,
	CellsToRemove const& remove,
	LinearGaussImproved *slopes,
	AMRExtensiveUpdater* eu) :
	refine_(refine), remove_(remove), scu_(vector<string>()),
	seu_(),interp_(slopes),
	eu_(eu)
{
	if (!eu)
		eu_ = &seu_;
}

void NonConservativeAMR::UpdateCellsRefine(Tessellation &tess,
	OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
	vector<Extensive> &extensives, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames,
	CacheData const& cd,
	PhysicalGeometry const& /*pg*/)const
{
	vector<size_t> ToRefine = refine_.ToRefine(tess, cells, time,tracerstickernames);
#ifndef RICH_MPI
	if (ToRefine.empty())
		return;
#endif // RICH_MPI
	sort(ToRefine.begin(), ToRefine.end());
	vector<std::pair<size_t, Vector2D> > NewPoints;
	vector<Vector2D> Moved;
#ifdef RICH_MPI
	vector<Vector2D> chull;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	ConvexHull(chull, proctess,rank);
#endif
	GetNewPoints(ToRefine, tess, NewPoints, Moved, obc
#ifdef RICH_MPI
		, chull
#endif
		);
	size_t N = static_cast<size_t>(tess.GetPointNo());
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);
	extensives.resize(N + NewPoints.size());

	for (size_t i = 0; i < NewPoints.size(); ++i)
	{
		cor.push_back(NewPoints[i].second);
		cells.push_back(cells[NewPoints[i].first]);
		if (interp_ != 0)
			interp_->GetSlopesUnlimited().push_back(interp_->GetSlopesUnlimited()[NewPoints[i].first]);
	}
	// Rebuild tessellation
#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif
	cd.reset();
	// Recalcualte extensives
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, cd.volumes[i],tracerstickernames);
#ifdef RICH_MPI
	MPI_exchange_data(tess, cells, true);
#endif
}

void NonConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames,
	CacheData const& cd,
	PhysicalGeometry const& /*pg*/)const
{
	std::pair<vector<size_t>, vector<double> > ToRemovepair = remove_.ToRemove(tess, cells, time,tracerstickernames);
	vector<size_t> ToRemove = ToRemovepair.first;
#ifndef RICH_MPI
	if (ToRemove.empty())
		return;
#endif // RICH_MPI
	size_t N = static_cast<size_t>(tess.GetPointNo());
	// Rebuild tessellation
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);

	RemoveVector(cor, ToRemove);
	RemoveVector(cells, ToRemove);

#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif
	cd.reset();
	// Recalcualte extensives
	extensives.resize(cells.size());
	for (size_t i = 0; i < extensives.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, cd.volumes[i],tracerstickernames);
#ifdef RICH_MPI
	MPI_exchange_data(tess, cells, true);
#endif
}

void NonConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames(), sim.getCacheData(), sim.getPhysicalGeometry());
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames(), sim.getCacheData(), sim.getPhysicalGeometry());
	// redo cache data
	sim.getCacheData().reset();
}

#ifdef RICH_MPI
vector<size_t> AMR::RemoveNearBoundaryPoints(vector<size_t> const&ToRemove,
	Tessellation const& tess, vector<double> &merits)const
{
	vector<size_t> res;
	int N = tess.GetPointNo();
	vector<double> merittemp;
	for (size_t i = 0; i < ToRemove.size(); ++i)
	{
		bool good = true;
		vector<int> neigh = tess.GetNeighbors(static_cast<int>(ToRemove[i]));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			if (tess.GetOriginalIndex(neigh[j]) >= N)
			{
				good = false;
				break;
			}
			else
			{
				if (tess.GetOriginalIndex(neigh[j]) == static_cast<int>(ToRemove[i]))
					continue;
				vector<int> neigh2 = tess.GetNeighbors(neigh[j]);
				for (size_t k = 0; k < neigh2.size(); ++k)
				{
					if (tess.GetOriginalIndex(neigh2[k]) == neigh[j])
						continue;
					else
					{
						if (tess.GetOriginalIndex(neigh2[k]) >= N)
						{
							good = false;
							break;
						}
					}
				}
			}
		}
		if (good)
		{
			res.push_back(ToRemove[i]);
			if(!merits.empty())
				merittemp.push_back(merits.at(i));
		}
	}
	merits = merittemp;
	return res;
}
#endif
