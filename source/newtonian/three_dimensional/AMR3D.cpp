#include "AMR3D.hpp"
#include "../../3D/GeometryCommon/Voronoi3D.hpp"
#include "../../misc/utils.hpp"

namespace
{
	void CleanOuterPoints(vector<size_t> &neigh, vector<size_t> &nneigh, Tessellation3D const&tess)
	{
		nneigh=RemoveList(nneigh, neigh);
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

	void RemoveBadAspectRatio(Tessellation3D const& tess, vector<double> &merits, vector<size_t> &toremove)
	{
		vector<size_t> remove_res;
		vector<double> merit_res;
		size_t N = toremove.size();
		for (size_t i = 0; i < N; ++i)
		{
			double v = tess.GetVolume(toremove[i]);
			vector<size_t> const& faces = tess.GetCellFaces(toremove[i]);
			size_t Nfaces = faces.size();
			double A = 0;
			for (size_t j = 0; j < Nfaces; ++j)
			{
				A += tess.GetArea(faces[j]);
			}			
			if ((A*A*A) < (v*v * 300))
			{
				remove_res.push_back(toremove[i]);
				merit_res.push_back(merits[i]);
			}
		}
		merits = merit_res;
		toremove = remove_res;
	}

#ifdef RICH_MPI
	void RemoveNeighbors()
	{
		bal
	}

#endif

	vector<size_t> RemoveNeighbors(vector<double> const& merits, vector<size_t> const& candidates, 
		Tessellation3D const& tess)
	{
		vector<size_t> result;
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
				result.push_back(candidates[i]);
			}
		}
		return result;
	}

	Vector3D GetNewPoint(Tessellation3D const& tess, vector<size_t> const& neigh, size_t index)
	{
		size_t Nneigh = neigh.size();
		assert(Nneigh > 0);
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
		return point*0.85 + 0.15*tess.GetMeshPoint(neigh[max_loc]);
	}

	/*void GetNeighborPoints(Tessellation3D const& tess, size_t index, vector<size_t> &neigh, vector<size_t> &nneigh)
	{
		tess.GetNeighbors(index, neigh);
		tess.GetNeighborNeighbors(nneigh, index);
	}*/

	void BuildLocalVoronoi(Tessellation3D &local, Tessellation3D const& tess, vector<size_t> const& real_points, 
		vector<size_t> const& nneigh,Vector3D const& newpoint,size_t torefine)
	{
		size_t Nreal = real_points.size();
		vector<Vector3D> points,ghost;
		points.reserve(Nreal + nneigh.size() + 1);
		points.push_back(newpoint);
		points.push_back(tess.GetMeshPoint(torefine));
		for (size_t i = 0; i < Nreal; ++i)
			points.push_back(tess.GetMeshPoint(real_points[i]));
		Nreal = nneigh.size();
		for (size_t i = 0; i < Nreal; ++i)
			if(nneigh[i]!=torefine)
				ghost.push_back(tess.GetMeshPoint(nneigh[i]));
		size_t counter = 2;
		local.BuildNoBox(points, ghost, 0);
	}

	bool IsRefineConserving(Tessellation3D const& local,size_t Nreal)
	{
		vector<size_t> const& faces = local.GetCellFaces(0);
		size_t Nfaces = faces.size();
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (local.GetFaceNeighbors(faces[i]).second >= Nreal && !local.IsPointOutsideBox(local.GetFaceNeighbors(faces[i]).second))
				return false;
		}
		return true;
	}

	Conserved3D CalcNewExtensives(Tessellation3D const& tess, Tessellation3D const& local, size_t torefine, vector<size_t> const& neigh,
		vector<ComputationalCell3D> const& cells,EquationOfState const& eos,TracerStickerNames const& tsn,
		vector<Conserved3D> &extensives)
	{
		Conserved3D res;
		PrimitiveToConserved(cells[torefine], tess.GetVolume(torefine),res,eos,tsn);
		Conserved3D newpoint(res);
		PrimitiveToConserved(cells[torefine], local.GetVolume(1), extensives[torefine], eos, tsn);
		newpoint -= extensives[torefine];
		assert(newpoint.mass > 0);
		size_t Nreal = neigh.size();
		for (size_t i = 0; i < Nreal; ++i)
		{
			if (tess.GetPointNo() <= neigh[i])
				continue;
			double dV = tess.GetVolume(neigh[i]) - local.GetVolume(i + 2);
			assert(dV > -1e-10*tess.GetVolume(neigh[i]));
			PrimitiveToConserved(cells[neigh[i]], dV, res, eos, tsn);
			PrimitiveToConserved(cells[neigh[i]], local.GetVolume(i + 2), extensives[neigh[i]], eos, tsn);
			newpoint += res;
		}
		return newpoint;
	}

	void GetBadIndeces(Tessellation3D const& tess, vector<size_t> const& neigh, size_t torefine, vector<size_t> &bad_faces)
	{
		bad_faces.clear();
		vector<size_t> const& faces = tess.GetCellFaces(torefine);
		size_t Nfaces = faces.size();
		for (size_t i = 0; i < Nfaces; ++i)
			bad_faces.push_back(faces[i]);
		size_t Nneigh = neigh.size();
		size_t Norg = tess.GetPointNo();
		for (size_t i = 0; i < Nneigh; ++i)
		{
			if (neigh[i] >= Norg)
				continue;
			vector<size_t> const& neigh_faces = tess.GetCellFaces(neigh[i]);
			size_t Nfaces2 = neigh_faces.size();
			for (size_t j = 0; j < Nfaces2; ++j)
			{
				std::pair<size_t, size_t> const& temp = tess.GetFaceNeighbors(neigh_faces[j]);
				size_t other = temp.first == neigh[i] ? temp.second : temp.first;
				if (std::binary_search(neigh.begin(), neigh.end(), other) || tess.IsPointOutsideBox(other))
					bad_faces.push_back(neigh_faces[j]);
			}
		}
		sort(bad_faces.begin(), bad_faces.end());
		bad_faces = unique(bad_faces);
	}

	void FixVoronoi(Tessellation3D &local, Tessellation3D &tess, vector<size_t> &neigh,size_t torefine,
		double &vol,Vector3D &CM,vector<size_t> const& nneigh,vector<Vector3D> &newboundary,vector<size_t> &newboundary_faces,
		size_t Nsplit)
	{
		// neigh is sorted
		size_t Norg = tess.GetPointNo();

		size_t Nlocal = neigh.size() + 2;
		size_t Nall = Nlocal + nneigh.size();
		vector<std::pair<size_t,size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		full_cellfaces.resize(full_cellfaces.size() + 1);
		vector<double> & full_area = tess.GetAllArea();
		size_t Nfaceslocal = localfaceneigh.size();

		vector<Vector3D>& full_vertices = tess.GetFacePoints();
		size_t Nvert = full_vertices.size();

		tess.GetAllCellFaces()[torefine].clear();
		vector<vector<size_t> > & all_cell_faces = tess.GetAllCellFaces();
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (tess.IsGhostPoint(neigh[i]))
				continue;
			vector<size_t> faces = tess.GetCellFaces(neigh[i]);
			vector<size_t> toremove;
			for (size_t j = 0; j < faces.size(); ++j)
			{
				std::pair<size_t, size_t> face_neigh = tess.GetFaceNeighbors(faces[j]);
				size_t other = face_neigh.first == neigh[i] ? face_neigh.second : face_neigh.first;
				if (tess.IsPointOutsideBox(other) || std::binary_search(neigh.begin(), neigh.end(), other) || other==torefine)
					toremove.push_back(j);
			}
			RemoveVector(all_cell_faces[neigh[i]], toremove);
		}

		vector<size_t> temp;
		size_t Ntotal = tess.GetMeshPoints().size();
		std::pair<size_t, size_t> new_face_neigh;
		for (size_t i = 0; i < Nfaceslocal; ++i)
		{
			if (localfaceneigh[i].second < Nlocal || ((localfaceneigh[i].first < Nlocal) && local.IsPointOutsideBox(localfaceneigh[i].second)))
			{
				if (localfaceneigh[i].first == 0)
					new_face_neigh.first = Ntotal;
				else
					if (localfaceneigh[i].first == 1)
						new_face_neigh.first = torefine;
					else
						new_face_neigh.first = neigh[localfaceneigh[i].first - 2];				
				if (localfaceneigh[i].second == 1)
					new_face_neigh.second = torefine;
				else
					if (localfaceneigh[i].second < Nlocal)
					{
						new_face_neigh.second = neigh[localfaceneigh[i].second - 2];
						if (tess.IsGhostPoint(new_face_neigh.first) && tess.IsGhostPoint(new_face_neigh.second))
							continue;
					}
					else
						if (localfaceneigh[i].second < Nall+4)
						{
							new_face_neigh.second = nneigh[localfaceneigh[i].second - Nlocal-4];
							if (tess.IsGhostPoint(new_face_neigh.first) && tess.IsGhostPoint(new_face_neigh.second))
								continue;
						}
						else
						{
							new_face_neigh.second = Ntotal + newboundary.size();
							if (tess.IsGhostPoint(new_face_neigh.first) && tess.IsGhostPoint(new_face_neigh.second))
								continue;
							newboundary.push_back(local.GetMeshPoint(localfaceneigh[i].second));
							newboundary_faces.push_back(full_faceneigh.size());
						}
				temp = local.GetPointsInFace(i);
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
				full_cellfaces.at(new_face_neigh.first).push_back(full_faceneigh.size() - 1);
				if(new_face_neigh.second<Norg)
					full_cellfaces.at(new_face_neigh.second).push_back(full_faceneigh.size() - 1);
				else
					if (new_face_neigh.second == Ntotal)
					{
						full_cellfaces.resize(Norg + Nsplit);
						full_cellfaces.back().push_back(full_faceneigh.size() - 1);
					}
				full_area.push_back(local.GetArea(i));
			}
		}
		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (neigh[i] < Norg)
			{
				tess.GetAllVolumes()[neigh[i]] = local.GetVolume(i+2);
				tess.GetAllCM()[neigh[i]] = local.GetCellCM(i + 2);
			}
		}
		tess.GetAllVolumes()[torefine] = local.GetVolume(1);
		vol = local.GetVolume(0);
		CM = local.GetCellCM(0);
		tess.GetAllCM()[torefine] = local.GetCellCM(1);
		full_vertices.insert(full_vertices.end(),local.GetFacePoints().begin(), local.GetFacePoints().end());


		for (size_t i = 0; i < neigh.size(); ++i)
		{
			if (tess.GetPointNo() <= neigh[i])
				continue;

			vector<size_t> t, t2;
			tess.GetNeighbors(neigh[i], t);
			local.GetNeighbors(i + 2, t2);
			assert(t.size() == t2.size());
		}
	}

	void FixBadIndeces(Tessellation3D &tess, vector<size_t> const& bad_indeces,size_t Norg2,size_t new_add)
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
		for (size_t i = Norg2; i < Norg; ++i)
		{
			vector<size_t> const& faces = tess.GetAllCellFaces()[i];
			size_t nfaces = faces.size();
			for (size_t j = 0; j < nfaces; ++j)
			{
				std::pair<size_t, size_t> neigh = tess.GetAllFaceNeighbors()[faces[j]];
				neigh.second -= new_add;
			}
		}
	}
}

AMRCellUpdater3D::~AMRCellUpdater3D(void) {}

AMRExtensiveUpdater3D::~AMRExtensiveUpdater3D(void){}

Conserved3D SimpleAMRExtensiveUpdater3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
	double volume, TracerStickerNames const& tracerstickernames) const
{
	Conserved3D res;
	const double mass = volume*cell.density;
	res.mass = mass;
	res.energy = eos.dp2e(cell.density, cell.pressure, cell.tracers)*mass +
		0.5*mass*ScalarProd(cell.velocity, cell.velocity);
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
	res.pressure = eos.de2p(res.density, extensive.energy / extensive.mass - 0.5*ScalarProd(res.velocity, res.velocity));
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

AMR3D::AMR3D(EquationOfState const& eos,CellsToRefine3D const& refine, CellsToRemove3D const& remove, LinearGauss3D *slopes, AMRCellUpdater3D* cu,
	AMRExtensiveUpdater3D* eu) :eos_(eos),refine_(refine), remove_(remove), interp_(slopes), cu_(cu), eu_(eu)
{
	if (!cu)
		cu_ = &scu_;
	if (!eu)
		eu_ = &seu_;
}

void AMR3D::UpdateCellsRefine(Tessellation3D &tess, vector<ComputationalCell3D> &cells, EquationOfState const& eos,
	vector<Conserved3D> &extensives, double time,
#ifdef RICH_MPI
	Tessellation3D const& proctess,
#endif
	TracerStickerNames const& tracerstickernames)const
{
	size_t Norg = tess.GetPointNo();
	std::pair<vector<size_t>,vector<double> > ToRefine = refine_.ToRefine(tess, cells, time, tracerstickernames);
	vector<size_t> indeces;
	sort_index(ToRefine.first, indeces);
	sort(ToRefine.first.begin(), ToRefine.first.end());
	ToRefine.second=VectorValues(ToRefine.second, indeces);
	RemoveBadAspectRatio(tess, ToRefine.second, ToRefine.first);
	ToRefine.first = RemoveNeighbors(ToRefine.second, ToRefine.first, tess);
	if (ToRefine.first.empty())
		return;
	size_t Nsplit = ToRefine.first.size();
	Voronoi3D vlocal(tess.GetBoxCoordinates().first,tess.GetBoxCoordinates().second);
	vector<size_t> neigh, nneigh,bad_faces,all_bad_faces,refined, newboundary_faces;
	double newvol;
	Vector3D newCM;
	vector<Vector3D> newpoints, newCMs,newboundary;
	vector<double> newvols;
	vector<Conserved3D> newcells;
	size_t new_add = tess.GetMeshPoints().size() - Norg;
	for (size_t i = 0; i < Nsplit; ++i)
	{
		tess.GetNeighbors(ToRefine.first[i], neigh);
		sort(neigh.begin(), neigh.end());
		tess.GetNeighborNeighbors(nneigh, ToRefine.first[i]);
		RemoveVal(nneigh, ToRefine.first[i]);
		//CleanOuterPoints(neigh, nneigh, tess);
		Vector3D NewPoint = GetNewPoint(tess, neigh, ToRefine.first[i]);
		BuildLocalVoronoi(vlocal, tess, neigh, nneigh, NewPoint, ToRefine.first[i]);
		if (IsRefineConserving(vlocal, neigh.size() + 2))
		{
			Conserved3D newcell = CalcNewExtensives(tess, vlocal, ToRefine.first[i], neigh,
				cells, eos, tracerstickernames, extensives);
			refined.push_back(ToRefine.first[i]);
			GetBadIndeces(tess, neigh, ToRefine.first[i], bad_faces);
			FixVoronoi(vlocal, tess, neigh, ToRefine.first[i], newvol, newCM,nneigh,newboundary,newboundary_faces,
				newpoints.size()+1);

			tess.GetMeshPoints().push_back(NewPoint);
			newvols.push_back(newvol);
			newCMs.push_back(newCM);
			newpoints.push_back(NewPoint);
			newcells.push_back(newcell);
			all_bad_faces.insert(all_bad_faces.end(), bad_faces.begin(), bad_faces.end());
		}
	}
	sort(all_bad_faces.begin(), all_bad_faces.end());
	all_bad_faces = unique(all_bad_faces);
	// Fix the old data and insert new data
	extensives.insert(extensives.begin() + Norg, newcells.begin(), newcells.end());
	vector<Vector3D> & allpoints = tess.GetMeshPoints();
	size_t Ntemp = allpoints.size();
	allpoints.insert(allpoints.begin() + Norg, newpoints.begin(), newpoints.end());
	allpoints.resize(Ntemp);
	allpoints.insert(allpoints.end(), newboundary.begin(), newboundary.end());
	for (size_t i = 0; i < newboundary_faces.size(); ++i)
		tess.GetAllFaceNeighbors()[newboundary_faces[i]].second = Ntemp + i;
	size_t &Norg2 = tess.GetPointNo();
	Norg2 += newpoints.size();
	vector<double> & allvol = tess.GetAllVolumes();
	allvol.insert(allvol.begin() + Norg, newvols.begin(), newvols.end());
	vector<Vector3D> & allCM = tess.GetAllCM();
	allCM.insert(allCM.begin() + Norg, newCMs.begin(), newCMs.end());
	cells.resize(Norg);
	// do new cells
	for (size_t i = 0; i < newCMs.size(); ++i)
		cells.push_back(cu_->ConvertExtensiveToPrimitve3D(extensives[i + Norg], eos, newvols[i], cells[refined[i]], tracerstickernames));
	RemoveVector(tess.GetAllFaceNeighbors(), all_bad_faces);
	RemoveVector(tess.GetAllArea(),all_bad_faces);
	RemoveVector(tess.GetAllPointsInFace(), all_bad_faces);
	// Fix all the indeces
	FixBadIndeces(tess, all_bad_faces,Norg,new_add);
}

void AMR3D::operator() (HDSim3D &sim)
{
	UpdateCellsRefine(sim.getTesselation(), sim.getCells(), eos_, sim.getExtensives(), sim.GetTime(), sim.GetTracerStickerNames());
}

AMR3D::~AMR3D(void){}