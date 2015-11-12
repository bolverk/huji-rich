#include "amr.hpp"
#include "../../tessellation/VoronoiMesh.hpp"

AMRCellUpdater::~AMRCellUpdater(){}

AMRExtensiveUpdater::~AMRExtensiveUpdater(){}

CellsToRemove::~CellsToRemove(){}

CellsToRefine::~CellsToRefine(){}

AMR::~AMR(void){}

Extensive SimpleAMRExtensiveUpdater::ConvertPrimitveToExtensive(const ComputationalCell& cell, const EquationOfState& eos,
	double volume) const
{
	Extensive res;
	const double mass = volume*cell.density;
	res.mass = mass;
	res.energy = eos.dp2e(cell.density, cell.pressure, cell.tracers)*mass +
		0.5*mass*ScalarProd(cell.velocity, cell.velocity);
	res.momentum = mass*cell.velocity;
	for (boost::container::flat_map<std::string, double>::const_iterator it =
		cell.tracers.begin();
		it != cell.tracers.end(); ++it)
		res.tracers[it->first] = (it->second)*mass;
	return res;
}

ComputationalCell SimpleAMRCellUpdater::ConvertExtensiveToPrimitve(const Extensive& extensive, const EquationOfState& eos,
	double volume, ComputationalCell const& old_cell) const
{
	ComputationalCell res;
	const double vol_inv = 1.0 / volume;
	res.density = extensive.mass*vol_inv;
	res.velocity = extensive.momentum/extensive.mass;
	res.pressure = eos.de2p(res.density, extensive.energy / extensive.mass - 0.5*ScalarProd(res.velocity, res.velocity));
	for (boost::container::flat_map<std::string, double>::const_iterator it =
		extensive.tracers.begin();
		it != extensive.tracers.end(); ++it)
		res.tracers[it->first] = (it->second)/extensive.mass;
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
		for (size_t i = 0; i < poly0.size(); ++i)
			maxi = std::max(maxi, std::max(std::abs(poly0[i].x), std::abs(poly0[i].y)));
		for (size_t i = 0; i < poly1.size(); ++i)
			maxi = std::max(maxi, std::max(std::abs(poly1[i].x), std::abs(poly1[i].y)));
		int maxscale = static_cast<int>(log10(maxi) + 9);

		subj[0].resize(poly0.size());
		clip[0].resize(poly1.size());
		for (size_t i = 0; i < poly0.size(); ++i)
			subj[0][i] = IntPoint(static_cast<cInt>(poly0[i].x*pow(10.0, 18 - maxscale)), static_cast<cInt>(poly0[i].y*pow(10.0, 18 - maxscale)));
		for (size_t i = 0; i < poly1.size(); ++i)
			clip[0][i] = IntPoint(static_cast<cInt>(poly1[i].x*pow(10.0, 18 - maxscale)), static_cast<cInt>(poly1[i].y*pow(10.0, 18 - maxscale)));

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

void AMR::GetNewPoints(vector<size_t> const& ToRefine, Tessellation const& tess,
	vector<std::pair<size_t, Vector2D> > &NewPoints, vector<Vector2D> &Moved,
	OuterBoundary const& obc
#ifdef RICH_MPI
	,vector<Vector2D> const& proc_chull
#endif
	)const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	NewPoints.clear();
	NewPoints.reserve(ToRefine.size() * 7);
	Moved.clear();
	Moved.reserve(ToRefine.size() * 7);
	for (size_t i = 0; i < ToRefine.size(); ++i)
	{
		vector<int> neigh = tess.GetNeighbors(static_cast<int>(ToRefine[i]));
		Vector2D const& mypoint = tess.GetMeshPoint(static_cast<int>(ToRefine[i]));
		const double R = tess.GetWidth(static_cast<int>(ToRefine[i]));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			Vector2D const& otherpoint = tess.GetMeshPoint(neigh[j]);
			if (otherpoint.distance(mypoint) < 1.75*R)
				continue;
			Vector2D candidate = 0.75*mypoint + 0.25*otherpoint;
#ifdef RICH_MPI
			if (!PointInCell(proc_chull, candidate))
				continue;
#endif
			if (static_cast<size_t>(neigh[j]) < N)
			{
				NewPoints.push_back(std::pair<size_t, Vector2D>(ToRefine[i],candidate));
				Moved.push_back(Vector2D());
			}
			else
			{
				Moved.push_back(FixInDomain(obc, candidate));
				NewPoints.push_back(std::pair<size_t, Vector2D>(ToRefine[i], candidate));
			}
		}
	}
}

ConservativeAMR::ConservativeAMR
(CellsToRefine const& refine,
	CellsToRemove const& remove,
	LinearGaussImproved *slopes,
	AMRCellUpdater* cu,
	AMRExtensiveUpdater* eu) :
	refine_(refine),
	remove_(remove),
	scu_(),
	seu_(),
	cu_(cu),
	eu_(eu),
	interp_(slopes)
{
  if(!cu){
    assert(!eu);
    cu_ = &scu_;
    eu_ = &seu_;
  }
}

void ConservativeAMR::UpdateCellsRefine
(Tessellation &tess,
 OuterBoundary const& obc,
 vector<ComputationalCell> &cells,
 EquationOfState const& eos,
 vector<Extensive> &extensives,
 double time
#ifdef RICH_MPI
	, Tessellation const& proctess
#endif
	)const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	// Find the primitive of each point and the new location
	vector<std::pair<size_t, Vector2D> > NewPoints;
	vector<Vector2D> Moved;
	vector<size_t> ToRefine = refine_.ToRefine(tess, cells, time);
#ifndef RICH_MPI
	if (ToRefine.empty())
		return;
#endif // RICH_MPI
#ifdef RICH_MPI
	ToRefine = RemoveNearBoundaryPoints(ToRefine, tess);
	vector<Vector2D> chull;
	const boost::mpi::communicator world;
	ConvexHull(chull, proctess, world.rank());
#endif
	GetNewPoints(ToRefine, tess, NewPoints, Moved, obc
#ifdef RICH_MPI
		,chull
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
#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif

	size_t location = 0;
	for (size_t i = 0; i < ToRefine.size(); ++i)
	{
		vector<int> neigh = oldtess->GetNeighbors(static_cast<int>(ToRefine[i]));
		vector<int> real_neigh;
		vector<vector<Vector2D> > Chull;
		vector<Vector2D> temp;
		ConvexHull(temp, *oldtess, static_cast<int>(ToRefine[i]));
		Chull.push_back(temp);
		real_neigh.push_back(static_cast<int>(ToRefine[i]));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			if (neigh[j] < static_cast<int>(N))
			{
				real_neigh.push_back(neigh[j]);
				ConvexHull(temp, *oldtess, neigh[j]);
				Chull.push_back(temp);
			}
			else
			{
				// Is it not a rigid wall?
				if (oldtess->GetOriginalIndex(neigh[j]) != static_cast<int>(ToRefine[i]))
				{
					real_neigh.push_back(oldtess->GetOriginalIndex(neigh[j]));
					ConvexHull(temp, *oldtess, real_neigh.front());
					temp = temp + (oldtess->GetMeshPoint(neigh[j]) -
						oldtess->GetMeshPoint(oldtess->GetOriginalIndex(neigh[j])));
					Chull.push_back(temp);
				}
			}
		}
		while (location<NewPoints.size()&&NewPoints[location].first == ToRefine[i])
		{
			double TotalVolume = 0;
			Extensive NewExtensive(extensives[0].tracers);
			for (size_t j = 0; j < Chull.size(); ++j)
			{
				ConvexHull(temp, tess, static_cast<int>(N + location));
				double v = AreaOverlap(temp, Chull[j]);
				NewExtensive += eu_->ConvertPrimitveToExtensive(cells[static_cast<size_t>(real_neigh[j])], eos, v);
				TotalVolume += v;
			}
			cells.push_back(cu_->ConvertExtensiveToPrimitve(NewExtensive, eos, TotalVolume, cells[ToRefine[i]]));
			if (interp_ != 0)
				interp_->GetSlopesUnlimited().push_back(interp_->GetSlopesUnlimited()[ToRefine[i]]);
			++location;
		}
	}
	extensives.resize(N + NewPoints.size());
	for (size_t i = 0; i < N + NewPoints.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)));
}

void ConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells,vector<Extensive> &extensives,
	EquationOfState const& eos,double time
#ifdef RICH_MPI
	, Tessellation const& proctess
#endif
	)const
{
	size_t N = static_cast<size_t>(tess.GetPointNo());
	// Rebuild tessellation
	vector<Vector2D> cor = tess.GetMeshPoints();
	cor.resize(N);
	cells.resize(N);
	std::pair<vector<size_t>,vector<double> > ToRemovepair = remove_.ToRemove(tess, cells, time);
	vector<size_t> ToRemove = ToRemovepair.first;
#ifndef RICH_MPI
	if (ToRemove.empty())
		return;
#endif // RICH_MPI
	ToRemove = RemoveNeighbors(ToRemovepair.second, ToRemovepair.first, tess);
#ifdef RICH_MPI
	ToRemove = RemoveNearBoundaryPoints(ToRemove, tess);
#endif
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
		Extensive oldcell = extensives[ToRemove[i]];
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
				extensives[static_cast<size_t>(oldtess->GetOriginalIndex(neigh[j]))] += (v/TotalV)*extensives[ToRemove[i]];
			}
		}
	}
	RemoveVector(extensives,ToRemove);
	RemoveVector(cells, ToRemove);
	N = static_cast<size_t>(tess.GetPointNo());
	for (size_t i = 0; i < N; ++i)
		cells[i] = cu_->ConvertExtensiveToPrimitve(extensives[i], eos, tess.GetVolume(static_cast<int>(i)), cells[i]);
}

NonConservativeAMR::NonConservativeAMR
  (CellsToRefine const& refine,
   CellsToRemove const& remove,
   AMRExtensiveUpdater* eu):
    refine_(refine), remove_(remove), scu_(),
	seu_(),
	eu_(eu)
{
  if(!eu)
    eu_ = &seu_;
}

void NonConservativeAMR::UpdateCellsRefine(Tessellation &tess,
	OuterBoundary const& obc, vector<ComputationalCell> &cells, EquationOfState const& eos,
	vector<Extensive> &extensives,double time
#ifdef RICH_MPI
	, Tessellation const& proctess
#endif
	)const
{
	vector<size_t> ToRefine = refine_.ToRefine(tess, cells, time);
#ifndef RICH_MPI
	if (ToRefine.empty())
		return;
#endif // RICH_MPI
	vector<std::pair<size_t, Vector2D> > NewPoints;
	vector<Vector2D> Moved;
#ifdef RICH_MPI
	vector<Vector2D> chull;
	const boost::mpi::communicator world;
	ConvexHull(chull, proctess, world.rank());
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
	extensives.resize(N+ToRefine.size());

	for (size_t i = 0; i < ToRefine.size(); ++i)
	{
		cor.push_back(NewPoints[i].second);
		cells.push_back(cells[NewPoints[i].first]);
	}
	// Rebuild tessellation
#ifdef RICH_MPI
	tess.Update(cor, proctess);
#else
	tess.Update(cor);
#endif

	// Recalcualte extensives
	for (size_t i = 0; i < N + ToRefine.size(); ++i)
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)));
}

void NonConservativeAMR::UpdateCellsRemove(Tessellation &tess,
	OuterBoundary const& /*obc*/, vector<ComputationalCell> &cells, vector<Extensive> &extensives,
	EquationOfState const& eos,double time
#ifdef RICH_MPI
	, Tessellation const& proctess
#endif
	)const
{
	std::pair<vector<size_t>,vector<double> > ToRemovepair = remove_.ToRemove(tess, cells, time);
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
		extensives[i] = eu_->ConvertPrimitveToExtensive(cells[i], eos, tess.GetVolume(static_cast<int>(i)));
}

void ConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(),sim.getTime()
#ifdef RICH_MPI
		, sim.GetProcTessellation()
#endif
		);
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(),sim.getTime()
#ifdef RICH_MPI
		, sim.GetProcTessellation()
#endif
		);
	// redo cache data
	sim.getCacheData().reset();
}

void NonConservativeAMR::operator()(hdsim &sim)
{
	UpdateCellsRefine(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getEos(),
		sim.getAllExtensives(), sim.getTime()
#ifdef RICH_MPI
		, sim.GetProcTessellation()
#endif
		);
	UpdateCellsRemove(sim.getTessellation(), sim.getOuterBoundary(), sim.getAllCells(), sim.getAllExtensives(),
		sim.getEos(), sim.getTime()
#ifdef RICH_MPI
		, sim.GetProcTessellation()
#endif
		);
	// redo cache data
	sim.getCacheData().reset();
}

vector<size_t> ConservativeAMR::RemoveNeighbors
(vector<double> const& merits, vector<size_t> const&
candidates, Tessellation const& tess) const
{
	vector<size_t> result;
	if (merits.size() != candidates.size())
		throw UniversalError("Merits and Candidates don't have same size in RemoveNeighbors");
	// Make sure there are no neighbors
	vector<size_t> bad_neigh;
	for (size_t i = 0; i<merits.size(); ++i)
	{
		bool good = true;
		vector<int> neigh = tess.GetNeighbors(static_cast<int>(candidates[i]));
		size_t nneigh = neigh.size();
		if (!good)
			continue;
		if (find(bad_neigh.begin(), bad_neigh.end(), candidates[i]) !=
			bad_neigh.end())
			good = false;
		else
		{
			for (size_t j = 0; j<nneigh; ++j)
			{
				if (binary_search(candidates.begin(), candidates.end(), neigh[j]))
				{
				  if (merits[i]<merits[static_cast<size_t>(lower_bound(candidates.begin(),candidates.end(), neigh[j]) - candidates.begin())])
					{
						good = false;
						break;
					}
				  if (fabs(merits[i] - merits[static_cast<size_t>(lower_bound(candidates.begin(),candidates.end(), neigh[j]) 
										  - candidates.begin())])<1e-9)
					{
						if (find(bad_neigh.begin(), bad_neigh.end(), neigh[j]) == bad_neigh.end())
							bad_neigh.push_back(static_cast<size_t>(neigh[j]));
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

#ifdef RICH_MPI
vector<size_t> AMR::RemoveNearBoundaryPoints(vector<size_t> const& candidates, Tessellation const& tess)const
{
	vector<size_t> res;
	int N = tess.GetPointNo();
	for (size_t i = 0; i < candidates.size(); ++i)
	{
		bool good = true;
		vector<int> neigh = tess.GetNeighbors(static_cast<int>(candidates[i]));
		for (size_t j = 0; j < neigh.size(); ++j)
		{
			if (tess.GetOriginalIndex(neigh[j]) >= N)
			{
				good = false;
				break;
			}
		}
		if (good)
			res.push_back(candidates[i]);
	}
	return res;
}
#endif
