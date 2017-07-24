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
		assert(Nbad > 1);
		vector<size_t> new_indeces(Nbad),bad_index(Nbad),total_remove;
		for (size_t i = 0; i < Nbad; ++i)
			bad_index[i] = i;
		while (bad_index.size() > 1)
		{
			vector<size_t> toremove;
			for (size_t i = 0; i < bad_index.size(); ++i)
			{
				if (abs(points[bad_indeces[bad_index[0]]] - points[bad_indeces[bad_index[i]]]) < 0.00001*R)
				{
					new_indeces[bad_index[i]] = bad_indeces[bad_index[0]];
					toremove.push_back(i);
					if (i > 0)
						total_remove.push_back(bad_indeces[bad_index[i]]);
				}
			}
			RemoveVector(bad_index, toremove);
			assert(bad_index.size() != 1);
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
