#include "Intersection3D.hpp"
#include <algorithm>
#include "../../misc/utils.hpp"

namespace
{/*
	class my_cmp
	{ 
	private:
		const Vector3D point_;

	public:
		my_cmp(Vector3D const& p) :point_(p) {}
		
		bool operator()(Vector3D const& v1, Vector3D const& v2)
		{
			double R = 1e-6*std::max(abs(v1), abs(v2));
			if (v1.x < (v2.x - R))
				return true;
			if (v1.x > (v2.x + R))
				return false;
			if (v1.y < (v2.y - R))
				return true;
			if (v1.y > (v2.y + R))
				return false;
			if (v1.z < (v2.z - R))
				return true;
			else
				return false;
		}
	};

*/
	void CleanDuplicates(vector<size_t> &all_indeces,vector<vector<int> > &face_inds,double R,vector<Vector3D> const& points)
	{
		size_t Nindeces = all_indeces.size();
		size_t Nfaces = face_inds.size();
		vector<size_t> counter(Nindeces, 0);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t Ninface = face_inds[i].size();
			for (size_t j = 0; j < Ninface; ++j)
			{
				vector<size_t>::const_iterator it = binary_find(all_indeces.begin(), all_indeces.end(),
					static_cast<size_t>(face_inds[i][j]));
				assert(it != all_indeces.end());
				++counter[static_cast<size_t>(it - all_indeces.begin())];
			}
		}
		vector<size_t> bad_indeces;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] < 3)
			{
				bad_indeces.push_back(all_indeces[i]);
			}
		}
		size_t Nbad = bad_indeces.size();
		if (Nbad < 1)
			return;

		vector<size_t> new_indeces(Nbad),bad_index(Nbad),total_remove;
		for (size_t i = 0; i < Nbad; ++i)
			bad_index[i] = i;
		while (bad_index.size() > 0)
		{
			vector<size_t> toremove;
			for (size_t i = 0; i < bad_index.size(); ++i)
			{
				if (abs(points[bad_indeces[bad_index[0]]] - points[bad_indeces[bad_index[i]]]) < 0.0003*R)
				{
					new_indeces[bad_index[i]] = bad_indeces[bad_index[0]];
					toremove.push_back(i);
					if (i > 0)
						total_remove.push_back(bad_indeces[bad_index[i]]);
				}
			}
			RemoveVector(bad_index, toremove);
		}
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t Ninface = face_inds[i].size();
			for (size_t j = 0; j < Ninface; ++j)
			{
				vector<size_t>::const_iterator it = binary_find(bad_indeces.begin(), bad_indeces.end(),
					static_cast<size_t>(face_inds[i][j]));
				if (it != bad_indeces.end())
				{
					face_inds[i][j] = static_cast<int>(new_indeces[static_cast<size_t>(it - bad_indeces.begin())]);
				}
			}
		}
		std::sort(total_remove.begin(), total_remove.end());
		all_indeces = RemoveList(all_indeces, total_remove);
	}

	void CleanDuplicates2(vector<size_t> &all_indeces, vector<vector<int> > &face_inds, double R, 
		vector<Vector3D> const& points,size_t index,Tessellation3D const& tess)
	{
		size_t Nindeces = all_indeces.size();
		size_t Nfaces = face_inds.size();
		vector<size_t> counter(Nindeces, 0);
		vector<size_t> const& faces = tess.GetCellFaces(index);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t Ninface = face_inds[i].size();
			for (size_t j = 0; j < Ninface; ++j)
			{
				vector<size_t>::const_iterator it = binary_find(all_indeces.begin(), all_indeces.end(),
					static_cast<size_t>(face_inds[i][j]));
				assert(it != all_indeces.end());
				++counter[static_cast<size_t>(it - all_indeces.begin())];
			}
		}
		vector<size_t> bad_indeces;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] < 3)
			{
				bad_indeces.push_back(all_indeces[i]);
				assert(counter[i] == 2);
			}
		}
		size_t Nbad = bad_indeces.size();
		if (Nbad < 1)
			return;
		// Find relevant faces
		vector<size_t> rel_faces(Nbad);
		for (size_t j = 0; j < Nbad; ++j)
		{
			size_t min_loc = 0;
			double min_dist = R * 1000;
			for (size_t i = 0; i < Nfaces; ++i)
			{
				if (std::find(face_inds[i].begin(), face_inds[i].end(), bad_indeces[j]) == face_inds[i].end())
				{
					size_t other = tess.GetFaceNeighbors(faces[i]).first == index ?
						tess.GetFaceNeighbors(faces[i]).second : tess.GetFaceNeighbors(faces[i]).first;
					double rnew = abs(tess.GetMeshPoint(other) - points[bad_indeces[j]]);
					if (rnew < min_dist)
					{
						min_dist = rnew;
						min_loc = i;
					}
				}
			}
			rel_faces[j] = min_loc;
		}
		// Add the vertex to the face
		for (size_t j = 0; j < Nbad; ++j)
		{
			Vector3D line = points[bad_indeces[j]] - tess.FaceCM(faces[rel_faces[j]]);
			line *= 1.0 / abs(line);
			size_t Ninface = face_inds[rel_faces[j]].size();
			size_t min_loc = 0;
			double max_angle = -2.0;
			for (size_t i = 0; i < Ninface; ++i)
			{
				Vector3D line2 = points[face_inds[rel_faces[j]][i]] - tess.FaceCM(faces[rel_faces[j]]);
				line2 *= 1.0 / abs(line2);
				double tempangle = ScalarProd(line2, line);
				if (tempangle > max_angle)
				{
					max_angle = tempangle;
					min_loc = i;
				}
			}
			Vector3D c_product1 = CrossProduct(line, points[face_inds[rel_faces[j]][min_loc]] 
				- tess.FaceCM(faces[rel_faces[j]]));
			Vector3D c_product2 = CrossProduct(line, points[face_inds[rel_faces[j]][(min_loc+1)%Ninface]]
				- tess.FaceCM(faces[rel_faces[j]]));
			if (ScalarProd(c_product1, c_product2) > 0)
				face_inds[rel_faces[j]].insert(face_inds[rel_faces[j]].begin() + (min_loc + Ninface - 1) % Ninface,
					static_cast<int>(bad_indeces[j]));
			else
				face_inds[rel_faces[j]].insert(face_inds[rel_faces[j]].begin() + (min_loc + 1) % Ninface,
					static_cast<int>(bad_indeces[j]));
		}
	}
}

void GetPlanes(vector<r3d_plane> &res, Tessellation3D const& tess, size_t index)
{
	res.clear();
	r3d_plane plane;
	vector<size_t> const& newfaces = tess.GetCellFaces(index);
	size_t nfaces = newfaces.size();
	res.resize(nfaces);
	for (size_t i = 0; i < nfaces; ++i)
	{
		std::pair<size_t, size_t> neigh = tess.GetFaceNeighbors(newfaces[i]);
		Vector3D norm = normalize(tess.GetMeshPoint(neigh.first) - tess.GetMeshPoint(neigh.second));
		if (neigh.second == index)
			norm *= -1;
		plane.n.xyz[0] = norm.x;
		plane.n.xyz[1] = norm.y;
		plane.n.xyz[2] = norm.z;
		Vector3D const& face_cm = tess.FaceCM(newfaces[i]);
		plane.d = -1 * (norm.x*face_cm.x + norm.y*face_cm.y + norm.z*face_cm.z);
		res[i] = plane;
	}
}

std::pair<bool, double> PolyhedraIntersection(Tessellation3D const & oldtess, Tessellation3D const & newtess, 
	size_t newcell, size_t oldcell,r3d_poly &poly,vector<Vector3D> &vtemp,vector<size_t> &itemp,vector<size_t> &all_indeces,
	vector<vector<int> > &faceinds, vector<r3d_plane> *planes)
{
	vtemp.clear();
	itemp.clear();
	all_indeces.clear();
	faceinds.clear();
	vector<size_t> const& oldfaces = oldtess.GetCellFaces(oldcell);
	size_t nfaces = oldfaces.size();
	faceinds.resize(nfaces);
	for (size_t i = 0; i < nfaces; ++i)
	{
		itemp = oldtess.GetPointsInFace(oldfaces[i]);
		faceinds[i].resize(itemp.size());
		for (size_t j = 0; j < itemp.size(); ++j)
			faceinds[i][j] = static_cast<int>(itemp[j]);
		all_indeces.insert(all_indeces.end(), itemp.begin(), itemp.end());
	}
	vector<int> numvertsperface(nfaces);
	std::sort(all_indeces.begin(), all_indeces.end());
	all_indeces = unique(all_indeces);
	vector<Vector3D> const& all_vertices = oldtess.GetFacePoints();
	// make sure no duplicate points
	CleanDuplicates(all_indeces, faceinds, oldtess.GetWidth(oldcell), all_vertices);
	CleanDuplicates2(all_indeces, faceinds, oldtess.GetWidth(oldcell), all_vertices, oldcell, oldtess);
	for (size_t i = 0; i < nfaces; ++i)
	{
		numvertsperface[i] = static_cast<int>(faceinds[i].size());
		for (size_t j = 0; j < faceinds[i].size(); ++j)
		{
			size_t index = static_cast<size_t>(binary_find(all_indeces.begin(), all_indeces.end(), static_cast<size_t>(
				faceinds[i][j])) - all_indeces.begin());
			faceinds[i][j] = static_cast<int>(index);
		}
		if (oldtess.GetFaceNeighbors(oldfaces[i]).second == oldcell)
			FlipVector(faceinds[i]);
	}
	size_t npoints = all_indeces.size();
	vector<r3d_rvec3> points(npoints);
	r3d_rvec3 point;
	for (size_t i = 0; i < npoints; ++i)
	{
		point.xyz[0] = all_vertices[all_indeces[i]].x;
		point.xyz[1] = all_vertices[all_indeces[i]].y;
		point.xyz[2] = all_vertices[all_indeces[i]].z;
		points[i] = point;
	}
	vector<int*> ptrs(faceinds.size());
	for (size_t i = 0, e = ptrs.size(); i<e; ++i) 
	{
		ptrs[i] = &(faceinds[i][0]);
	}
	r3d_init_poly(&poly, &points[0], static_cast<r3d_int>(points.size()), &ptrs[0], &numvertsperface[0], 
		static_cast<r3d_int>(numvertsperface.size()));
	int test = r3d_is_good(&poly);
	assert(test == 1);
	bool allocated = false;
	if (planes == 0)
	{
		allocated = true;
		planes = new vector<r3d_plane>;
		GetPlanes(*planes, newtess, newcell);
	}
	r3d_clip(&poly, &(planes->at(0)), static_cast<int>(planes->size()));
	std::pair<bool, double> res(false,0);
	if (poly.nverts == 0)
		res.first = false;
	else
	{
		res.first = true;
		// calc volume of intersection
		double m[1];
		m[0] = 0;
		r3d_reduce(&poly, m, 0);
		res.second = std::abs(m[0]);
	}
	if (allocated)
		delete planes;
	return res;
}
