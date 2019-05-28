#include "amr.hpp"
#include "../../tessellation/VoronoiMesh.hpp"

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
	double volume, TracerStickerNames const& /*tracerstickernames*/) const
{
	Extensive res;
	const double mass = volume*cell.density;
	res.mass = mass;
	res.energy = eos.dp2e(cell.density, cell.pressure, cell.tracers)*mass +
		0.5*mass*ScalarProd(cell.velocity, cell.velocity);
	res.momentum = mass*cell.velocity;
	size_t N = cell.tracers.size();
//	res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell.tracers[i] * mass;
	return res;
}

SimpleAMRCellUpdater::SimpleAMRCellUpdater(vector<string> toskip) :toskip_(toskip) {}

ComputationalCell SimpleAMRCellUpdater::ConvertExtensiveToPrimitve(const Extensive& extensive, const EquationOfState& eos,
	double volume, ComputationalCell const& old_cell,TracerStickerNames const& tracerstickernames) const
{
	for (size_t i = 0; i < toskip_.size(); ++i)
		if(*safe_retrieve(old_cell.stickers.begin(), tracerstickernames.sticker_names.begin(),
			tracerstickernames.sticker_names.end(),toskip_[i]))
			return old_cell;
	ComputationalCell res;
	const double vol_inv = 1.0 / volume;
	res.density = extensive.mass*vol_inv;
	res.velocity = extensive.momentum / extensive.mass;
	res.pressure = eos.de2p(res.density, extensive.energy / extensive.mass - 0.5*ScalarProd(res.velocity, res.velocity));
	size_t N = extensive.tracers.size();
//	res.tracers.resize(N);
	for (size_t i = 0; i < N;++i)
		res.tracers[i]=extensive.tracers[i] / extensive.mass;
	res.stickers = old_cell.stickers;
	return res;
}

namespace
{
	double AreaOverlap(vector<Vector2D> const& poly0, vector<Vector2D> const& poly1)
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
		PolygonOverlap polyoverlap;
		if (!solution.empty())
		{
			if (solution[0].size() > 2)
			{
				vector<Vector2D> inter(solution[0].size());
				double factor = pow(10.0, maxscale - 18);
				for (size_t i = 0; i < solution[0].size(); ++i)
					inter[i].Set(static_cast<double>(solution[0][i].X) * factor, static_cast<double>(solution[0][i].Y) * factor);
				return polyoverlap.PolyArea(inter);
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

	vector<size_t> RemoveNeighbors(vector<double> const& merits, vector<size_t> const&
		candidates, Tessellation const& tess)
	{
		vector<size_t> result;
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
				result.push_back(candidates[i]);
			}
		}
		return result;
	}

	Extensive GetNewExtensive(vector<Extensive> const& extensives,Tessellation const& tess,size_t N,size_t location,
		vector<Vector2D> const& moved,vector<int> const& real_neigh,vector<vector<Vector2D> >  const& Chull,
		vector<ComputationalCell> const& cells,EquationOfState const& eos, AMRExtensiveUpdater const& eu,
		double &TotalVolume,Tessellation const& oldtess,bool periodic,TracerStickerNames const& tracerstickernames)
	{
		Extensive NewExtensive(extensives[0].tracers);
		vector<Vector2D> temp;
		ConvexHull(temp, tess, static_cast<int>(N + location));
		temp = temp + moved[location];
		const double vv = tess.GetVolume(static_cast<int>(N + location));
		stack<int> tocheck;
		stack<vector<Vector2D> > tocheck_hull;
		vector<int> checked(real_neigh);
		vector<int> neightemp;
		for (size_t j = 0; j < Chull.size(); ++j)
		{
			double v = AreaOverlap(temp, Chull[j]);
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
			double v = AreaOverlap(temp, chull);
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
		TracerStickerNames const& tracerstickernames)
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
				cells, eos, eu,TotalVolume,oldtess,periodic,tracerstickernames);
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
	TracerStickerNames const& tracerstickernames)const
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

#ifdef debug_amr
	WriteTess(tess, "c:/vnew.h5");
#endif


	size_t location = 0;
	for (size_t i = 0; i < ToRefine.size(); ++i)
		ConservedSingleCell(*oldtess, tess, ToRefine[i], location, interp_, eos, cells, NewPoints, extensives, *eu_,
			*cu_,Moved,periodic_,tracerstickernames);
	extensives.resize(N + NewPoints.size());
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)),tracerstickernames);
}

void ConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames)const
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
#ifdef RICH_MPI
	ToRemovepair.first = RemoveNearBoundaryPoints(ToRemovepair.first, tess,ToRemovepair.second);
#endif

	vector<size_t> ToRemove = RemoveNeighbors(ToRemovepair.second, ToRemovepair.first, tess);
	// save copy of old tessellation
	boost::scoped_ptr<Tessellation> oldtess(tess.clone());

	RemoveVector(cor, ToRemove);

#ifdef debug_amr
	WriteTess(tess, "c:/vold.h5");
#endif

#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif

#ifdef debug_amr
	WriteTess(tess, "c:/vnew.h5");
#endif

	// Fix the extensives
	vector<Vector2D> temp;
	for (size_t i = 0; i < ToRemove.size(); ++i)
	{
		vector<int> neigh = oldtess->GetNeighbors(static_cast<int>(ToRemove[i]));
		vector<Vector2D> chull;
		ConvexHull(chull, *oldtess, static_cast<int>(ToRemove[i]));
		const double TotalV = oldtess->GetVolume(static_cast<int>(ToRemove[i]));
		temp.clear();
		double dv = 0;
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			size_t toadd = static_cast<size_t>(lower_bound(ToRemove.begin(), ToRemove.end(),
				oldtess->GetOriginalIndex(neigh[j])) - ToRemove.begin());
			temp.clear();
			if (neigh[j] < static_cast<int>(N))
			{
				ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(neigh[j]) - toadd));
			}
			else
			{
				if (oldtess->GetOriginalIndex(neigh[j]) != static_cast<int>(ToRemove[i]))
				{
					ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(
						oldtess->GetOriginalIndex(neigh[j]) )- toadd));
					temp = temp + (oldtess->GetMeshPoint(neigh[j]) -
						oldtess->GetMeshPoint(oldtess->GetOriginalIndex(neigh[j])));
				}
			}
			if (!temp.empty())
			{
				const double v = AreaOverlap(chull, temp);
				dv += v;
				extensives[static_cast<size_t>(oldtess->GetOriginalIndex(neigh[j]))] += (v / TotalV)*extensives[ToRemove[i]];
			}
		}
		if (dv > (1 + 1e-5)*TotalV || dv < (1 - 1e-5)*TotalV)
		{
			std::cout << "Real volume: " << TotalV << " AMR volume: " << dv << std::endl;
			throw UniversalError("Not same volume in amr remove");
		}
	}
	RemoveVector(extensives, ToRemove);
	RemoveVector(cells, ToRemove);

	N = static_cast<size_t>(tess.GetPointNo());
	for (size_t i = 0; i < N; ++i)
		cells[i] = cu_->ConvertExtensiveToPrimitve(extensives[i], eos, tess.GetVolume(static_cast<int>(i)), cells[i],tracerstickernames);
}

void ConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
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
	TracerStickerNames const& tracerstickernames)const
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

	// Recalcualte extensives
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)),tracerstickernames);
}

void NonConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames)const
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
	// Recalcualte extensives
	extensives.resize(cells.size());
	for (size_t i = 0; i < extensives.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)),tracerstickernames);
}

void NonConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
	// redo cache data
	sim.getCacheData().reset();
}

ConservativeAMROld::ConservativeAMROld
(CellsToRefine const& refine,
	CellsToRemove const& remove,
	LinearGaussImproved *slopes,
	AMRCellUpdater* cu,
	AMRExtensiveUpdater* eu) :
	refine_(refine),
	remove_(remove),
	scu_(vector<string>()),
	seu_(),
	cu_(cu),
	eu_(eu),
	interp_(slopes)
{
	if (!cu)
		cu_ = &scu_;
	if(!eu)
		eu_ = &seu_;
}

void ConservativeAMROld::UpdateCellsRefine
(Tessellation &tess,
	OuterBoundary const& obc,
	vector<ComputationalCell> &cells,
	EquationOfState const& eos,
	vector<Extensive> &extensives,
	double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames)const
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
#ifdef RICH_MPI
	vector<double> merit;
	ToRefine = RemoveNearBoundaryPoints(ToRefine, tess,merit);
#endif
	GetNewPoints2(ToRefine, tess, NewPoints, Moved, obc);

	// Rebuild tessellation
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);
	if (interp_ != 0)
		interp_->GetSlopesUnlimited().resize(N);

	for (size_t i = 0; i < NewPoints.size(); ++i)
	{
		cor.push_back(NewPoints[i].second);
		cells.push_back(cells[NewPoints[i].first]);
		interp_->GetSlopesUnlimited().push_back(interp_->GetSlopesUnlimited().at(NewPoints[i].first));
	}
#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif

	extensives.resize(N + NewPoints.size());
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)),tracerstickernames);
}

void ConservativeAMROld::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos, double time,
#ifdef RICH_MPI
	Tessellation const& proctess,
#endif
	TracerStickerNames const& tracerstickernames)const
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
#ifdef RICH_MPI
	ToRemovepair.first = RemoveNearBoundaryPoints(ToRemovepair.first, tess,ToRemovepair.second);
#endif
	vector<size_t> ToRemove = RemoveNeighbors(ToRemovepair.second, ToRemovepair.first, tess);
	// save copy of old tessellation
	boost::scoped_ptr<Tessellation> oldtess(tess.clone());

	RemoveVector(cor, ToRemove);
#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif

	// Fix the extensives
	vector<Vector2D> temp;
	for (size_t i = 0; i < ToRemove.size(); ++i)
	{
		vector<int> neigh = oldtess->GetNeighbors(static_cast<int>(ToRemove[i]));
		vector<Vector2D> chull;
		ConvexHull(chull, *oldtess, static_cast<int>(ToRemove[i]));
		const double TotalV = oldtess->GetVolume(static_cast<int>(ToRemove[i]));
		//		Extensive oldcell = extensives[ToRemove[i]];
		temp.clear();
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			size_t toadd = static_cast<size_t>(lower_bound(ToRemove.begin(), ToRemove.end(), neigh[j]) - ToRemove.begin());
			if (neigh[j] < static_cast<int>(N))
			{
				ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(neigh[j]) - toadd));
			}
			else
			{
				if (oldtess->GetOriginalIndex(neigh[j]) != static_cast<int>(ToRemove[i]))
				{
					ConvexHull(temp, tess, static_cast<int>(static_cast<size_t>(neigh[j]) - toadd));
					temp = temp + (oldtess->GetMeshPoint(neigh[j]) -
						oldtess->GetMeshPoint(oldtess->GetOriginalIndex(neigh[j])));
				}
			}
			if (!temp.empty())
			{
				const double v = AreaOverlap(chull, temp);
				extensives[static_cast<size_t>(oldtess->GetOriginalIndex(neigh[j]))] += (v / TotalV)*extensives[ToRemove[i]];
			}
		}
	}
	RemoveVector(extensives, ToRemove);
	RemoveVector(cells, ToRemove);
	N = static_cast<size_t>(tess.GetPointNo());
	for (size_t i = 0; i < N; ++i)
		cells[i] = cu_->ConvertExtensiveToPrimitve(extensives[i], eos, tess.GetVolume(static_cast<int>(i)), cells[i],tracerstickernames);
}

void ConservativeAMROld::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime(),
#ifdef RICH_MPI
		sim.GetProcTessellation(),
#endif
		sim.GetTracerStickerNames());
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
