#include "Intersection3D.hpp"
#include <algorithm>
#include "../../misc/utils.hpp"
#include <array>

namespace
{
	double AbsDiff(r3d_rvec3 const& v1, r3d_rvec3 const& v2)
	{
		double res = 0;
		res += (v1.xyz[0] - v2.xyz[0])*(v1.xyz[0] - v2.xyz[0]);
		res += (v1.xyz[1] - v2.xyz[1])*(v1.xyz[1] - v2.xyz[1]);
		res += (v1.xyz[2] - v2.xyz[2])*(v1.xyz[2] - v2.xyz[2]);
		return std::sqrt(res);
	}

	void FixBadOrder(vector<vector<int> > & faceinds, vector<r3d_rvec3> &all_vert, double R)
	{
		size_t Nfaces = faceinds.size();
		size_t Npoints = all_vert.size();
		vector<vector<int> > right_neighbor(Npoints);
		vector<int > toremove;
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t NinFace = faceinds[i].size();
			for (size_t j = 0; j < NinFace; ++j)
			{
				if (right_neighbor[static_cast<size_t>(faceinds[i][j])].empty() ||
					std::find(right_neighbor[static_cast<size_t>(faceinds[i][j])].begin(),
						right_neighbor[static_cast<size_t>(faceinds[i][j])].end(), faceinds[i][(j + 1) % NinFace])
					== right_neighbor[static_cast<size_t>(faceinds[i][j])].end())
					right_neighbor[faceinds[i][j]].push_back(faceinds[i][(j + 1) % NinFace]);
				else
				{
					double Rdiff = AbsDiff(all_vert[faceinds[i][j]], all_vert[faceinds[i][(j + 1) % NinFace]]);
					if (Rdiff > 1e-3*R)
					{
						std::cout << "R " << R << " Rdiff " << Rdiff << " face " << i << " points " << all_vert[faceinds[i][j]].xyz[0]
							<< " " << all_vert[faceinds[i][j]].xyz[1] << " " << all_vert[faceinds[i][j]].xyz[2]
							<< " " << all_vert[faceinds[i][(j + 1) % NinFace]].xyz[0] << " " << all_vert[faceinds[i][(j + 1) % NinFace]].xyz[1]
							<< " " << all_vert[faceinds[i][(j + 1) % NinFace]].xyz[2] << std::endl;
					}
					toremove.push_back(faceinds[i][j]);
				}
			}
		}
		std::sort(toremove.begin(), toremove.end());
		toremove = unique(toremove);
		// remove points that have bad order (these should be duplicate points)
		for (size_t i = 0; i < Nfaces; ++i)
			faceinds[i] = RemoveList(faceinds[i], toremove);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			for (size_t j = 0; j < faceinds[i].size(); ++j)
			{
				int index_change = static_cast<int>(std::lower_bound(toremove.begin(), toremove.end(),
					faceinds[i][j]) - toremove.begin());
				faceinds[i][j] -= index_change;
			}
		}
		for (size_t i = 0; i < toremove.size(); ++i)
			RemoveVector(all_vert, toremove);
	}

	vector<size_t> GetBadPoints(vector<vector<int> > &face_inds, vector<size_t> &all_indeces)
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
		vector<size_t> bad_indeces, to_remove;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] < 3 && counter[i]>0)
				bad_indeces.push_back(all_indeces[i]);
			if (counter[i] == 0)
				to_remove.push_back(i);
		}
		vector<int> toremoveint;
		for (size_t i = 0; i < to_remove.size(); ++i)
			toremoveint.push_back(static_cast<int>(all_indeces[to_remove[i]]));
		std::sort(toremoveint.begin(), toremoveint.end());
		RemoveVector(all_indeces, to_remove);
		for (size_t i = 0; i < face_inds.size(); ++i)
			face_inds[i] = RemoveList(face_inds[i], toremoveint);
		std::sort(bad_indeces.begin(), bad_indeces.end());
		return bad_indeces;
	}

	void CleanDuplicates(vector<size_t> &all_indeces, vector<vector<int> > &face_inds, vector<Vector3D> const& points)
	{
		size_t Nindeces = all_indeces.size();
		vector<size_t> bad_indeces = GetBadPoints(face_inds, all_indeces);

		size_t Nbad = bad_indeces.size();
		if (Nbad < 1)
			return;

		vector<size_t> all_good = RemoveList(all_indeces, bad_indeces);
		all_good.insert(all_good.end(), bad_indeces.begin(), bad_indeces.end());

		vector<double> dist;
		size_t Nfaces = face_inds.size();
		for (size_t i = 0; i < Nbad; ++i)
		{
			size_t to_change = bad_indeces.back();
			bad_indeces.pop_back();
			all_good.pop_back();
			Nindeces = all_good.size();
			dist.resize(Nindeces);
			for (size_t j = 0; j < Nindeces; ++j)
				dist[j] = fastabs(points[all_good[j]] - points[to_change]);
			size_t min_loc = static_cast<size_t>(std::distance(dist.begin(), std::min_element(dist.begin(), dist.end())));
			size_t max_loc = static_cast<size_t>(std::distance(dist.begin(), std::max_element(dist.begin(), dist.end())));
			if (dist[min_loc] > dist[max_loc] * 1e-6)
				continue;
			for (size_t j = 0; j < Nfaces; ++j)
				for (size_t k = 0; k < face_inds[j].size(); ++k)
					if (face_inds[j][k] == static_cast<int>(to_change))
						face_inds[j][k] = static_cast<int>(all_good[min_loc]);
		}
		// clean all_indeces
		all_indeces.clear();
		for (size_t j = 0; j < Nfaces; ++j)
			for (size_t k = 0; k < face_inds[j].size(); ++k)
				all_indeces.push_back(static_cast<size_t>(face_inds[j][k]));
		std::sort(all_indeces.begin(), all_indeces.end());
		all_indeces = unique(all_indeces);
	}

	void RemoveTooMany(vector<vector<int> > &face_inds, vector<size_t> &too_many)
	{
		size_t Nface = face_inds.size();
		vector<size_t> to_remove;
		for (size_t k = 0; k < too_many.size(); ++k)
		{
			for (size_t i = 0; i < Nface; ++i)
			{
				to_remove.clear();
				size_t Npoints = face_inds[i].size();
				bool found = false;
				for (size_t j = 0; j < Npoints; ++j)
				{
					if (face_inds[i][j] == static_cast<int>(too_many[k]))
					{
						if (found)
							to_remove.push_back(j);
						found = true;
					}
				}
				if (!to_remove.empty())
					RemoveVector(face_inds[i], to_remove);
			}
		}
	}

	bool CleanDuplicates2(vector<size_t> &all_indeces, vector<vector<int> > &face_inds)
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
		vector<size_t> too_many;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] > 3)
				too_many.push_back(all_indeces[i]);
		}
		if (!too_many.empty())
			RemoveTooMany(face_inds, too_many);
		too_many.clear();
		for (size_t i = 0; i < Nfaces; ++i)
			if (face_inds[i].size() < 3)
				too_many.push_back(i);
		RemoveVector(face_inds, too_many);
		Nfaces = face_inds.size();
		counter.assign(Nindeces, 0);
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
		too_many.clear();
		bool res = true;
		for (size_t i = 0; i < Nindeces; ++i)
			if (counter[i] == 0)
				too_many.push_back(i);
			else
				if (counter[i] != 3)
					res = false;
		RemoveVector(all_indeces, too_many);
		return res;
	}

	void OrganizeToPolyData(vector<vector<int> > &faceinds, vector<int> &numvertsperface, vector<size_t> &all_indeces,
		vector<Vector3D> const& all_vertices, vector<r3d_rvec3> &points, vector<int*> &ptrs)
	{
		size_t nfaces = faceinds.size();
		for (size_t i = 0; i < nfaces; ++i)
		{
			numvertsperface[i] = static_cast<int>(faceinds[i].size());
			for (size_t j = 0; j < faceinds[i].size(); ++j)
			{
				size_t index = static_cast<size_t>(binary_find(all_indeces.begin(), all_indeces.end(), static_cast<size_t>(
					faceinds[i][j])) - all_indeces.begin());
				faceinds[i][j] = static_cast<int>(index);
			}
		}
		size_t npoints = all_indeces.size();
		r3d_rvec3 point;
		for (size_t i = 0; i < npoints; ++i)
		{
			point.xyz[0] = all_vertices[all_indeces[i]].x;
			point.xyz[1] = all_vertices[all_indeces[i]].y;
			point.xyz[2] = all_vertices[all_indeces[i]].z;
			points[i] = point;
		}
		for (size_t i = 0, e = ptrs.size(); i < e; ++i)
		{
			ptrs[i] = &(faceinds[i][0]);
		}
	}
}

void GetPlanes(vector<r3d_plane> &res, Tessellation3D const& tess, size_t index)
{
	res.clear();
	r3d_plane plane;
	face_vec const& newfaces = tess.GetCellFaces(index);
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

bool GetPoly(Tessellation3D const & oldtess, size_t oldcell, r3d_poly &poly, point_vec &itemp, vector<size_t> &all_indeces,
	vector<vector<int> > &faceinds)
{
	itemp.clear();
	all_indeces.clear();
	faceinds.clear();
	face_vec const& oldfaces = oldtess.GetCellFaces(oldcell);
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

	std::sort(all_indeces.begin(), all_indeces.end());
	all_indeces = unique(all_indeces);
	vector<Vector3D> const& all_vertices = oldtess.GetFacePoints();
	for (size_t i = 0; i < nfaces; ++i)
		if (oldtess.GetFaceNeighbors(oldfaces[i]).second == oldcell)
			FlipVector(faceinds[i]);

	// make sure no duplicate points
	CleanDuplicates(all_indeces, faceinds, all_vertices);
	CleanDuplicates2(all_indeces, faceinds);
	CleanDuplicates(all_indeces, faceinds, all_vertices);
	bool goodpoly = CleanDuplicates2(all_indeces, faceinds);
	if (!goodpoly)
	{
		std::cout << "Bad polygon in cell " << oldcell << std::endl;
		all_indeces.clear();
		for (size_t i = 0; i < nfaces; ++i)
		{
			itemp = oldtess.GetPointsInFace(oldfaces[i]);
			std::cout << "Face " << oldfaces[i] << " points:";
			for (size_t j = 0; j < itemp.size(); ++j)
				std::cout << " " << itemp[j] << " ";
			std::cout << std::endl;
			all_indeces.insert(all_indeces.end(), itemp.begin(), itemp.end());
		}
		sort(all_indeces.begin(), all_indeces.end());
		all_indeces = unique(all_indeces);
		for (size_t i = 0; i < all_indeces.size(); ++i)
		{
			Vector3D p = all_vertices[all_indeces[i]];
			std::cout << "Point " << all_indeces[i] << " " << p.x << " " << p.y << " " << p.z << std::endl;
		}
		return false;
	}
	nfaces = faceinds.size();
	size_t npoints = all_indeces.size();
	vector<int> numvertsperface(nfaces);
	vector<r3d_rvec3> points(npoints);
	vector<int*> ptrs(faceinds.size());
	OrganizeToPolyData(faceinds, numvertsperface, all_indeces, all_vertices, points, ptrs);
	if (points.size() > 256)
	{
		std::cout << "Too many points in cells refine " << points.size() << std::endl;
		return false;
	}
	r3d_init_poly(&poly, &points[0], static_cast<r3d_int>(points.size()), &ptrs[0], &numvertsperface[0],
		static_cast<r3d_int>(numvertsperface.size()));
	int test = r3d_is_good(&poly);
	if (test != 1)
	{
		FixBadOrder(faceinds, points, oldtess.GetWidth(oldcell));
		for (size_t i = 0; i < nfaces; ++i)
			numvertsperface[i] = static_cast<int>(faceinds[i].size());
		for (size_t i = 0, e = ptrs.size(); i < e; ++i)
			ptrs[i] = &(faceinds[i][0]);
		r3d_init_poly(&poly, &points[0], static_cast<r3d_int>(points.size()), &ptrs[0], &numvertsperface[0],
			static_cast<r3d_int>(numvertsperface.size()));
		test = r3d_is_good(&poly);
		if (test != 1)
		{
			std::cout << "Bad polygon in cell " << oldcell << std::endl;
			std::cout << "Given to poly build:" << std::endl;
			for (size_t i = 0; i < npoints; ++i)
			{
				std::cout << "Point " << all_indeces[i] << " " << all_vertices[all_indeces[i]].x <<
					" " << all_vertices[all_indeces[i]].y << " " << all_vertices[all_indeces[i]].z << std::endl;
			}
			for (size_t i = 0; i < ptrs.size(); ++i)
			{
				std::cout << "Face " << i;
				for (size_t j = 0; j < faceinds[i].size(); ++j)
					std::cout << " " << faceinds[i][j] << " ";
				std::cout << std::endl;
			}

			all_indeces.clear();
			for (size_t i = 0; i < nfaces; ++i)
			{
				itemp = oldtess.GetPointsInFace(oldfaces[i]);
				std::cout << "Face " << oldfaces[i] << " points:";
				for (size_t j = 0; j < itemp.size(); ++j)
					std::cout << " " << itemp[j] << " ";
				std::cout << std::endl;
				all_indeces.insert(all_indeces.end(), itemp.begin(), itemp.end());
			}
			sort(all_indeces.begin(), all_indeces.end());
			all_indeces = unique(all_indeces);
			for (size_t i = 0; i < all_indeces.size(); ++i)
			{
				Vector3D p = all_vertices[all_indeces[i]];
				std::cout << "Point " << all_indeces[i] << " " << p.x << " " << p.y << " " << p.z << std::endl;
			}
			return false;
		}
		return true;
	}
	return true;
}

std::pair<bool, std::array<double,4> > PolyhedraIntersection(Tessellation3D const & newtess, size_t newcell, r3d_poly &poly,
	vector<r3d_plane> *planes)
{
	bool allocated = false;
	if (planes == 0)
	{
		allocated = true;
		planes = new vector<r3d_plane>;
		GetPlanes(*planes, newtess, newcell);
	}
	r3d_clip(&poly, &(planes->at(0)), static_cast<int>(planes->size()));
	std::pair<bool, std::array<double, 4> > res;
	if (poly.nverts == 0)
		res.first = false;
	else
	{
		res.first = true;
		// calc volume of intersection and it's CM
		std::array<double, 4> m;
		m[0] = 0;
		r3d_reduce(&poly, &m[0], 1);
		m[0] = std::abs(m[0]);
		if (!(m[0] > 0))
			res.first = false;
		for (size_t j = 1; j < 4; ++j)
			m[j] /= m[0];
		res.second = m;
	}
	if (allocated)
		delete planes;
	return res;
}
