#include "AMR3D.hpp"
#include "../../3D/GeometryCommon/Voronoi3D.hpp"
#include "../../misc/utils.hpp"

//#define debug_amr 1

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

	void RemoveBadAspectRatio(Tessellation3D const& tess, vector<size_t> &toremove)
	{
		vector<size_t> remove_res;
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
			if (((A*A*A) < (v*v * 300))&&(abs(tess.GetMeshPoint(toremove[i])-tess.GetCellCM(toremove[i]))<
				0.2*tess.GetWidth(toremove[i])))
			{
				remove_res.push_back(toremove[i]);
			}
		}
		toremove = remove_res;
	}

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
				vector<Vector3D> const& vertices = tess.GetFacePoints();
				vector<size_t> indeces = tess.GetPointsInFace(faces[j]);
				Vector3D normal = CrossProduct(vertices[indeces[1]] - vertices[indeces[0]],
					vertices[indeces[2]] - vertices[indeces[0]]);
				normal *= (1.0 / abs(normal));
				if (tess.GetFaceNeighbors(faces[j]).second == i)
					normal *= -1;
				sum += normal*a;
			}
			vector<size_t> neigh;
			if (abs(sum) > 1e-5*A)
				tess.GetNeighbors(i, neigh);
			assert(abs(sum) < 1e-4*A);
			assert(PointInPoly(tess, tess.GetMeshPoint(i), i));
		}
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
		return point*0.999999 + (1- 0.999999)*tess.GetMeshPoint(neigh[max_loc]);
	}


	void BuildLocalVoronoi(Tessellation3D &local, Tessellation3D const& tess, vector<size_t> const& real_points, 
		Vector3D const& newpoint,size_t torefine)
	{
		size_t Nreal = real_points.size();
		vector<Vector3D> points,ghost;
		points.push_back(newpoint);
		points.push_back(tess.GetMeshPoint(torefine));
		ghost.reserve(Nreal);
		for (size_t i = 0; i < Nreal; ++i)
			ghost.push_back(tess.GetMeshPoint(real_points[i]));
		local.BuildNoBox(points, ghost, 0);
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

	void FixVoronoi(Tessellation3D &local, Tessellation3D &tess, vector<size_t> &neigh,size_t torefine,
		double &vol,Vector3D &CM,size_t Ntotal0,size_t index)
	{
		// neigh is sorted
		vector<std::pair<size_t, size_t> >const& localfaceneigh = local.GetAllFaceNeighbors();
		vector<std::pair<size_t, size_t> >& full_faceneigh = tess.GetAllFaceNeighbors();
		vector<vector<size_t> >& full_facepoints = tess.GetAllPointsInFace();
		vector<vector<size_t> >& full_cellfaces = tess.GetAllCellFaces();
		vector<double> & full_area = tess.GetAllArea();
		vector<Vector3D>& full_vertices = tess.GetFacePoints();

		size_t Norg = tess.GetPointNo();
		vector<size_t> faces = tess.GetCellFaces(torefine);
		size_t Nfaces = faces.size();
		// Remove old face reference
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t other = tess.GetFaceNeighbors(faces[i]).first==torefine ? tess.GetFaceNeighbors(faces[i]).second :
				tess.GetFaceNeighbors(faces[i]).first;
			if (tess.IsPointOutsideBox(other))
				continue;
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
		size_t Ntotal = tess.GetMeshPoints().size();
		size_t Nlocal = neigh.size() + 6;
		size_t Nvert = full_vertices.size();
		full_cellfaces.resize(Ntotal0 + index + 1);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			if (localfaceneigh[faces[i]].first == 0)
				new_face_neigh.first = Ntotal0+index;
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
				{
					new_face_neigh.second = Ntotal;
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
			if (new_face_neigh.first<Norg || ((new_face_neigh.first>=Ntotal0)&&(new_face_neigh.first<Ntotal0+index)))
				full_cellfaces.at(new_face_neigh.first).push_back(full_faceneigh.size() - 1);
			else
			{
				assert(new_face_neigh.first == Ntotal0+index);
				full_cellfaces.at(new_face_neigh.first).push_back(full_faceneigh.size() - 1);
			}

			if (new_face_neigh.second<Norg || ((new_face_neigh.second >= Ntotal0) && (new_face_neigh.second<=Ntotal0 + index)))
				full_cellfaces.at(new_face_neigh.second).push_back(full_faceneigh.size() - 1);
			full_area.push_back(local.GetArea(i));
		}
		full_vertices.insert(full_vertices.end(), local.GetFacePoints().begin(), local.GetFacePoints().end());
	}

	void FixBadIndeces(Tessellation3D &tess, vector<size_t> const& bad_indeces,size_t Nsplit,size_t Ntotal0)
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
					tess.GetAllFaceNeighbors()[i].first += Norg - Ntotal0-Nsplit;
			}
			else
			{
				if (tess.GetFaceNeighbors(i).first >= (Norg - Nsplit))
					tess.GetAllFaceNeighbors()[i].first += Nsplit;
			}
			if (tess.GetFaceNeighbors(i).second >= Ntotal0)
			{
				if (tess.GetFaceNeighbors(i).second < Ntotal0 + Nsplit)
					tess.GetAllFaceNeighbors()[i].second += Norg - Ntotal0-Nsplit;
			}
			else
			{
				if (tess.GetFaceNeighbors(i).second >= (Norg - Nsplit))
					tess.GetAllFaceNeighbors()[i].second += Nsplit;
			}

		}
	}
}

AMRCellUpdater3D::~AMRCellUpdater3D(void) {}

AMRExtensiveUpdater3D::~AMRExtensiveUpdater3D(void){}

Conserved3D SimpleAMRExtensiveUpdater3D::ConvertPrimitveToExtensive3D(const ComputationalCell3D& cell, const EquationOfState& eos,
	double volume, TracerStickerNames const& /*tracerstickernames*/) const
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
	Tessellation3D const& /*proctess*/,
#endif
	TracerStickerNames const& tracerstickernames)const
{
	size_t Norg = tess.GetPointNo();
	size_t Ntotal0 = tess.GetMeshPoints().size();
	extensives.resize(Norg);
	cells.resize(Norg);
	vector<size_t> ToRefine = refine_.ToRefine(tess, cells, time, tracerstickernames);
	vector<size_t> indeces;
	sort_index(ToRefine, indeces);
	sort(ToRefine.begin(), ToRefine.end());
	RemoveBadAspectRatio(tess, ToRefine);
	if (ToRefine.empty())
		return;
	size_t Nsplit = ToRefine.size();
	Voronoi3D vlocal(tess.GetBoxCoordinates().first,tess.GetBoxCoordinates().second);
	vector<size_t> neigh, bad_faces,all_bad_faces,refined, newboundary_faces;
	double newvol;
	Vector3D newCM;
	vector<Vector3D> newpoints, newCMs;
	vector<double> newvols;
	tess.GetMeshPoints().resize(Ntotal0 + Nsplit);
	for (size_t i = 0; i < Nsplit; ++i)
	{
		tess.GetNeighbors(ToRefine[i], neigh);
		sort(neigh.begin(), neigh.end());
		Vector3D NewPoint = GetNewPoint(tess, neigh, ToRefine[i]);
		BuildLocalVoronoi(vlocal, tess, neigh, NewPoint, ToRefine[i]);
		bad_faces = tess.GetCellFaces(ToRefine[i]);
		/////////////
		double oldv = tess.GetVolume(ToRefine[i]);
		double newv = vlocal.GetVolume(0) + vlocal.GetVolume(1);
		assert(oldv > 0.9999*newv&&newv > 0.9999*oldv);
		///////////
		FixVoronoi(vlocal, tess, neigh, ToRefine[i], newvol, newCM,Ntotal0,i);
		PrimitiveToConserved(cells[ToRefine[i]], tess.GetVolume(ToRefine[i]), extensives[ToRefine[i]], eos, tracerstickernames);
		tess.GetMeshPoints()[Ntotal0 + i] = NewPoint;
		newvols.push_back(newvol);
		newCMs.push_back(newCM);
		newpoints.push_back(NewPoint);
		all_bad_faces.insert(all_bad_faces.end(), bad_faces.begin(), bad_faces.end());
	}
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
		cells.push_back(cells[ToRefine[i]]);
		PrimitiveToConserved(cells[ToRefine[i]], tess.GetVolume(Norg+i), extensives[Norg+i], eos, tracerstickernames);
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
	vector<vector<size_t> > temp_cell_faces(tess.GetAllCellFaces().begin() + Ntotal0, tess.GetAllCellFaces().begin() + Nsplit +Ntotal0);
	std::copy(temp_cell_faces.begin(), temp_cell_faces.end(), tess.GetAllCellFaces().begin() + Norg);
	tess.GetAllCellFaces().resize(Norg2);

	RemoveVector(tess.GetAllFaceNeighbors(), all_bad_faces);
	RemoveVector(tess.GetAllArea(),all_bad_faces);
	RemoveVector(tess.GetAllPointsInFace(), all_bad_faces);
	// Fix all the indeces
	FixBadIndeces(tess, all_bad_faces,Nsplit,Ntotal0);
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
}

AMR3D::~AMR3D(void){}