#include "Intersection3D.hpp"
#include <algorithm>
#include "../../misc/utils.hpp"
#include <limits>
#include <boost/container/flat_map.hpp>

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

	void FixBadOrder(vector<vector<int> > & faceinds, vector<r3d_rvec3> &all_vert,double R)
	{
		size_t Nfaces = faceinds.size();
		size_t Npoints = all_vert.size();
		size_t smax = std::numeric_limits<size_t>::max();
		boost::container::flat_map<size_t,  size_t > right_neighbor;
		right_neighbor.reserve(Npoints);
		for (size_t i = 0; i < Npoints; ++i)
			right_neighbor.insert(std::pair<size_t,size_t>(i,smax));
		vector<int > toremove;
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t NinFace = faceinds[i].size();
			for (size_t j = 0; j < NinFace; ++j)
			{
				if (right_neighbor[faceinds[i][j]] == smax)
					right_neighbor[faceinds[i][j]] = faceinds[i][(j + 1) % NinFace];
				else
				{
					double Rdiff = AbsDiff(all_vert[faceinds[i][j]], all_vert[faceinds[i][(j + 1) % NinFace]]);
					if (Rdiff > 1e-3*R)
					{
						std::cout << "R " << R << " Rdiff " << Rdiff << std::endl;
					}
					assert(Rdiff < 1e-3*R);
					toremove.push_back(faceinds[i][j]);
				}
			}
		}
		std::sort(toremove.begin(), toremove.end());
		toremove = unique(toremove);
		// remove points that have bad order (these should be duplicate points)
		for (size_t i = 0; i < Nfaces; ++i)
			faceinds[i] = RemoveList(faceinds[i], toremove);
		for (size_t i = 0; i < toremove.size(); ++i)
			RemoveVector(all_vert, toremove);
	}

	void FixDegenerate(vector<vector<int> > & faceinds, size_t bad_point,vector<size_t> &all_indeces)
	{
		size_t Nfaces = faceinds.size();
		size_t bad_before = 0;
		size_t bad_after = 0;
		bool found = false;
		all_indeces.push_back(all_indeces.back() + 1);
		vector<size_t> containing_faces;
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t NinFace = faceinds[i].size();
			for (size_t j = 0; j < NinFace; ++j)
			{
				if (faceinds[i][j] == static_cast<int>(bad_point) && !found)
				{
					bad_after = faceinds[i][(j + 1) % NinFace];
					bad_before = faceinds[i][(j + NinFace - 1) % NinFace];
					faceinds[i][j] = static_cast<int>(all_indeces.back());
					found = true;
					break;
				}
				else
					if (faceinds[i][j] == static_cast<int>(bad_point))
					{
						containing_faces.push_back(i);
						break;
					}
			}
		}
		size_t counter = 0;
		for (size_t i = 0; i < containing_faces.size(); ++i)
		{
			size_t face = containing_faces[i];
			size_t NinFace = faceinds[face].size();
			for (size_t j = 0; j < NinFace; ++j)
			{
				if (faceinds[face][j] == static_cast<int>(bad_after))
				{
					faceinds[face].insert(faceinds[face].begin()+(j + 1) % NinFace, static_cast<int>(all_indeces.back()));
					++counter;
					break;
				}
				if (faceinds[face][j] == static_cast<int>(bad_before))
				{
					faceinds[face].insert(faceinds[face].begin()+j % NinFace, static_cast<int>(all_indeces.back()));
					++counter;
					break;
				}
			}
			if (counter == 2)
				break;
		}
	}

	vector<size_t> GetBadPoints(vector<vector<int> > const& face_inds,vector<size_t> &all_indeces)
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
		vector<size_t> bad_indeces,to_remove;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] < 3 &&counter[i]>0)
				bad_indeces.push_back(all_indeces[i]);
			if (counter[i] == 0)
				to_remove.push_back(i);
		}
		RemoveVector(all_indeces, to_remove);
		std::sort(bad_indeces.begin(), bad_indeces.end());
		return bad_indeces;
	}

	vector<vector<std::pair<size_t, size_t> > > GetBadFaces(vector<vector<int> > const& face_inds,vector<size_t> const&
		bad_indeces)
	{
		size_t Nbad = bad_indeces.size();
		size_t Nfaces = face_inds.size();
		vector<vector<std::pair<size_t, size_t> > > bad_faces(Nbad);
		for (size_t i = 0; i < Nfaces; ++i)
		{
			size_t Ninface = face_inds[i].size();
			for (size_t j = 0; j < Ninface; ++j)
			{
				vector<size_t>::const_iterator it = binary_find(bad_indeces.begin(), bad_indeces.end(),
					static_cast<size_t>(face_inds[i][j]));
				if (it != bad_indeces.end())
					bad_faces[static_cast<size_t>(it - bad_indeces.begin())].push_back(std::pair<size_t, size_t>(i, j));
			}
		}
		return bad_faces;
	}
	void CleanDuplicates(vector<size_t> &all_indeces,vector<vector<int> > &face_inds,double R,vector<Vector3D> const& points)
	{
		size_t Nindeces = all_indeces.size();
		vector<size_t> bad_indeces = GetBadPoints(face_inds, all_indeces);
		
		size_t Nbad = bad_indeces.size();
		if (Nbad < 1)
			return;
		

		vector<double> dist;
		vector<size_t> sort_ind;
		while (Nbad > 0)
		{
			Nindeces = all_indeces.size();
			dist.resize(Nindeces);
			vector<vector<std::pair<size_t, size_t> > > bad_faces = GetBadFaces(face_inds, bad_indeces);
			for (size_t i = 0; i < Nbad; ++i)
			{
				for (size_t j = 0; j < Nindeces; ++j)
				{
					if (bad_indeces[i] == all_indeces[j])
						dist[j] = R * 100;
					else
						dist[j] = abs(points[all_indeces[j]] - points[bad_indeces[i]]);
				}
				sort_index(dist, sort_ind);
				size_t value = std::min(bad_indeces[i], all_indeces[sort_ind[0]]);
				if (value < bad_indeces[i])
				{
					for (size_t j = 0; j < bad_faces[i].size(); ++j)
						face_inds[bad_faces[i][j].first][bad_faces[i][j].second] = static_cast<int>(value);
				}
				else
					if(!std::binary_search(bad_indeces.begin(),bad_indeces.end(),all_indeces[sort_ind[0]]))
						for (size_t j = 0; j < bad_faces[i].size(); ++j)
							face_inds[bad_faces[i][j].first][bad_faces[i][j].second] = static_cast<int>(all_indeces[sort_ind[0]]);
			}
			bad_indeces = GetBadPoints(face_inds, all_indeces);
			Nbad = bad_indeces.size();
		}
/*
		vector<size_t> new_indeces(Nbad),bad_index(Nbad),total_remove;
		for (size_t i = 0; i < Nbad; ++i)
			bad_index[i] = i;
		while (bad_index.size() > 0)
		{
			vector<size_t> toremove;
			for (size_t i = 0; i < bad_index.size(); ++i)
			{
				if (abs(points[bad_indeces[bad_index[0]]] - points[bad_indeces[bad_index[i]]]) < 0.0005*R)
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
		all_indeces = RemoveList(all_indeces, total_remove);*/
	}

	void RemoveTooMany(vector<vector<int> > &face_inds,vector<size_t> &too_many)
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

	void CleanDuplicates2(vector<size_t> &all_indeces, vector<vector<int> > &face_inds, vector<size_t> &degn_points)
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
		vector<size_t> bad_indeces,too_many;
		for (size_t i = 0; i < Nindeces; ++i)
		{
			if (counter[i] < 3)
			{
				bad_indeces.push_back(all_indeces[i]);
				//assert(counter[i] == 2);
			}
			if (counter[i] > 3)
			{
				too_many.push_back(all_indeces[i]);
			}
		}
		if (!too_many.empty())
			RemoveTooMany(face_inds, too_many);
	/*	size_t Nbad = bad_indeces.size();
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
				face_inds[rel_faces[j]].insert(face_inds[rel_faces[j]].begin() + min_loc,
					static_cast<int>(bad_indeces[j]));
			else
				face_inds[rel_faces[j]].insert(face_inds[rel_faces[j]].begin() + (min_loc + 1) % Ninface,
					static_cast<int>(bad_indeces[j]));
		}*/
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
		degn_points.clear();
		for (size_t i = 0; i < Nindeces; ++i)
		{
			assert(counter[i] > 2);
			if (counter[i] > 3) // degenerate case with a vertex shared by 4 faces, remove smallest face
			{
				FixDegenerate(face_inds, all_indeces[i], all_indeces);
				degn_points.push_back(all_indeces[i]);
			}
		}
	}

	void OrganizeToPolyData(vector<vector<int> > &faceinds, vector<int> &numvertsperface, vector<size_t> &all_indeces,
		vector<Vector3D> const& all_vertices, size_t Npoints, vector<size_t> const& degn_points,
		vector<r3d_rvec3> &points, vector<int*> &ptrs)
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
			if (i >= Npoints)
			{
				point.xyz[0] = all_vertices[degn_points[i - Npoints]].x;
				point.xyz[1] = all_vertices[degn_points[i - Npoints]].y;
				point.xyz[2] = all_vertices[degn_points[i - Npoints]].z;
			}
			else
			{
				point.xyz[0] = all_vertices[all_indeces[i]].x;
				point.xyz[1] = all_vertices[all_indeces[i]].y;
				point.xyz[2] = all_vertices[all_indeces[i]].z;
			}
			points[i] = point;
		}
		for (size_t i = 0, e = ptrs.size(); i<e; ++i)
		{
			ptrs[i] = &(faceinds[i][0]);
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

void GetPoly(Tessellation3D const & oldtess,size_t oldcell, r3d_poly &poly, vector<size_t> &itemp, vector<size_t> &all_indeces,
	vector<vector<int> > &faceinds)
{
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
	for (size_t i = 0; i < nfaces; ++i)
		if (oldtess.GetFaceNeighbors(oldfaces[i]).second == oldcell)
			FlipVector(faceinds[i]);
	// make sure no duplicate points
	CleanDuplicates(all_indeces, faceinds, oldtess.GetWidth(oldcell), all_vertices);
	size_t Npoints = all_indeces.size();
	vector<size_t> degn_points;
	CleanDuplicates2(all_indeces, faceinds, degn_points);
	size_t npoints = all_indeces.size();
	vector<r3d_rvec3> points(npoints);
	vector<int*> ptrs(faceinds.size());
	OrganizeToPolyData(faceinds, numvertsperface, all_indeces, all_vertices, Npoints, degn_points,
		points, ptrs);
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
				if (i >= Npoints)
					std::cout << "Point " << degn_points[i - Npoints] << " " << all_vertices[degn_points[i - Npoints]].x <<
					" " << all_vertices[degn_points[i - Npoints]].y << " " << all_vertices[degn_points[i - Npoints]].z << std::endl;
				else
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
		}
	}
	assert(test == 1);
}

std::pair<bool, double> PolyhedraIntersection(Tessellation3D const & newtess,size_t newcell,r3d_poly &poly,
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
