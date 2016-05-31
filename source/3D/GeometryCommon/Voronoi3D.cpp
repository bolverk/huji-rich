#include "Voronoi3D.hpp"
#include <algorithm>
#include <stack>
#include "Mat44.hpp"
#include "Intersections.hpp"
#include "utils.hpp"
#include <fstream>
#include <iostream>


namespace
{
	double CalcFaceArea(vector<size_t> const& indeces,vector<Vector3D> const& points)
	{
		size_t Nloop = indeces.size() - 2;
		Vector3D temp;
		for (size_t i = 0; i < Nloop; ++i)
			temp += CrossProduct(points[indeces[i + 2]] - points[indeces[0]], points[indeces[i + 1]] - points[indeces[0]]);
		return 0.5*abs(temp);
	}

	class VecCompare
	{
	public:
		double R;
		vector<Vector3D> const& points;
		VecCompare(vector<Vector3D> const& mesh) :points(mesh), R(0) {}
		bool operator() (size_t i, size_t j)
		{
			if (abs(points[i].x - points[j].x) > 1e-8*R)
				return points[i].x < points[j].x;
			else
				if (abs(points[i].y - points[j].y) > 1e-8*R)
					return points[i].y < points[j].y;
				else
					return points[i].z < points[j].z;
		}
	};


	void RemoveDuplicates(vector<Vector3D> const& points, vector<size_t> &indeces, double R,
		vector<size_t> &temp, VecCompare &compare)
	{
		compare.R = R;
		std::sort(indeces.begin(), indeces.end(), compare);
		double eps = 1e-8;
		size_t N = indeces.size();
		temp.clear();
		temp.push_back(indeces[0]);
		for (size_t i = 1; i < N; ++i)
			if (abs(points[indeces[i]] - points[temp.back()])>R*eps)
				temp.push_back(indeces[i]);
		indeces = temp;
	}

	double CalcRadius(boost::array<Vector3D, 4> const&points)
	{
#define v1 (points[0])
#define v2 (points[1])
#define v3 (points[2])
#define v4 (points[3])

		Mat44<double> m_a(v1.x, v1.y, v1.z, 1,
			v2.x, v2.y, v2.z, 1,
			v3.x, v3.y, v3.z, 1,
			v4.x, v4.y, v4.z, 1);
		double a = m_a.determinant();

		Mat44<double> m_Dx(ScalarProd(v1, v1), v1.y, v1.z, 1,
			ScalarProd(v2, v2), v2.y, v2.z, 1,
			ScalarProd(v3, v3), v3.y, v3.z, 1,
			ScalarProd(v4, v4), v4.y, v4.z, 1);
		double Dx = m_Dx.determinant();

		Mat44<double> m_Dy(ScalarProd(v1, v1), v1.x, v1.z, 1,
			ScalarProd(v2, v2), v2.x, v2.z, 1,
			ScalarProd(v3, v3), v3.x, v3.z, 1,
			ScalarProd(v4, v4), v4.x, v4.z, 1);
		double Dy = -m_Dy.determinant();

		Mat44<double> m_Dz(ScalarProd(v1, v1), v1.x, v1.y, 1,
			ScalarProd(v2, v2), v2.x, v2.y, 1,
			ScalarProd(v3, v3), v3.x, v3.y, 1,
			ScalarProd(v4, v4), v4.x, v4.y, 1);
		double Dz = m_Dz.determinant();

		Mat44<double> m_c(ScalarProd(v1, v1), v1.x, v1.y, v1.z,
			ScalarProd(v2, v2), v2.x, v2.y, v2.z,
			ScalarProd(v3, v3), v3.x, v3.y, v3.z,
			ScalarProd(v4, v4), v4.x, v4.y, v4.z);
		double c = m_c.determinant();

#undef v1
#undef v2
#undef v3
#undef v4
		return 0.5*sqrt(Dx*Dx + Dy*Dy + Dz*Dz - 4 * a*c) / std::abs(a);
	}


	void ConvexHull3D(vector<Vector3D> const& points, vector<size_t> &indeces, vector<size_t> & temp)
	{
		// Find center point in plane
		Vector3D center = (points[indeces[0]] + points[indeces[1]] + points[indeces[2]]);
		for (size_t i = 3; i < indeces.size(); ++i)
			center = center + points[indeces[i]];
		center = center / indeces.size();
		vector<double> angles(indeces.size() - 1);
		size_t Nloop = angles.size();
		Vector3D main_vector = points[indeces[Nloop]] - center;
		main_vector = main_vector / abs(main_vector);
		Vector3D normal = CrossProduct(main_vector, center - points[indeces[0]]);
		normal = normal / abs(normal);
		for (size_t i = 0; i < Nloop; ++i)
		{
			Vector3D other_vector = points[indeces[i]] - center;
			other_vector = other_vector / abs(other_vector);
			double temp = ScalarProd(main_vector, other_vector);
			if (ScalarProd(normal, CrossProduct(other_vector, main_vector)) < 0)
				temp += -2 - 2 * temp;
			angles[i] = temp;
		}
		temp = sort_index(angles);
		temp.push_back(Nloop);
		ReArrangeVector(indeces, temp);
	}

	vector<Face> BuildBox(Vector3D const& ll, Vector3D const& ur)
	{
		double dx = ur.x - ll.x;
		double dy = ur.y - ll.y;
		double dz = ur.z - ll.z;
		vector<Face> res(6);
		vector<Vector3D> points;
		points.push_back(ll);
		points.push_back(ll + Vector3D(dx, 0, 0));
		points.push_back(ll + Vector3D(dx, dy, 0));
		points.push_back(ll + Vector3D(0, dy, 0));
		points.push_back(ll + Vector3D(0, 0, dz));
		points.push_back(ll + Vector3D(dx, 0, dz));
		points.push_back(ll + Vector3D(dx, dy, dz));
		points.push_back(ll + Vector3D(0, dy, dz));
		points.push_back(ur);
		res[0].vertices.push_back(points[0]);
		res[0].vertices.push_back(points[1]);
		res[0].vertices.push_back(points[2]);
		res[0].vertices.push_back(points[3]);
		res[1].vertices.push_back(points[0]);
		res[1].vertices.push_back(points[4]);
		res[1].vertices.push_back(points[5]);
		res[1].vertices.push_back(points[1]);
		res[2].vertices.push_back(points[3]);
		res[2].vertices.push_back(points[7]);
		res[2].vertices.push_back(points[4]);
		res[2].vertices.push_back(points[0]);
		res[3].vertices.push_back(points[2]);
		res[3].vertices.push_back(points[6]);
		res[3].vertices.push_back(points[7]);
		res[3].vertices.push_back(points[3]);
		res[4].vertices.push_back(points[1]);
		res[4].vertices.push_back(points[5]);
		res[4].vertices.push_back(points[6]);
		res[4].vertices.push_back(points[2]);
		res[5].vertices.push_back(points[5]);
		res[5].vertices.push_back(points[4]);
		res[5].vertices.push_back(points[7]);
		res[5].vertices.push_back(points[6]);
		return res;
	}

	bool PointInVertices(b_array_4 const& points, size_t point)
	{
		return !(std::find(points.begin(), points.end(), point) == points.end());
	}

	Vector3D MirrorPoint(Face const& face, Vector3D const& point)
	{
		Vector3D normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
		normal = normal / abs(normal);
		return point - (2 * ScalarProd(point - face.vertices[0], normal))*normal;
	}

}

Voronoi3D::Voronoi3D()
{}

Voronoi3D::Voronoi3D(Vector3D const& ll, Vector3D const& ur) :ll_(ll), ur_(ur) {}

Tetrahedron::Tetrahedron(void) {}

void Voronoi3D::RunTetGen(vector<Vector3D> const& points, bool voronoi)
{
	tetin.deinitialize();
	tetout.deinitialize();
	tetin.firstnumber = 0;
	tetin.numberofpoints = points.size();
	size_t N = points.size();
	tetin.pointlist = new double[N * 3 + 12];
	for (size_t i = 0; i < N; ++i)
	{
		tetin.pointlist[i * 3] = points[i].x;
		tetin.pointlist[i * 3 + 1] = points[i].y;
		tetin.pointlist[i * 3 + 2] = points[i].z;
	}
	if (!voronoi)
	{
		Vector3D pmax(points[0]), pmin(points[0]);

		for (size_t i = 0; i < N; ++i)
		{
			pmax.x = std::max(pmax.x, points[i].x);
			pmax.y = std::max(pmax.y, points[i].y);
			pmax.z = std::max(pmax.z, points[i].z);
			pmin.x = std::min(pmin.x, points[i].x);
			pmin.y = std::min(pmin.y, points[i].y);
			pmin.z = std::min(pmin.z, points[i].z);
		}
		mesh_points_.reserve(points.size() + 4);
		mesh_points_ = points;
		// Build big tetrahedron
		mesh_points_.push_back(Vector3D(0.5*(pmax.x + pmin.x), 0.5*(pmax.y + pmin.y), pmax.z + 10 * (pmax.z - pmin.z)));
		mesh_points_.push_back(Vector3D(pmin.x - 10 * (pmax.x - pmin.x), pmin.y - 10 * (pmax.y - pmin.y),
			pmin.z - 100 * (pmax.z - pmin.z)));
		mesh_points_.push_back(Vector3D(pmax.x + 10 * (pmax.x - pmin.x), pmin.y - 10 * (pmax.y - pmin.y),
			pmin.z - 100 * (pmax.z - pmin.z)));
		mesh_points_.push_back(Vector3D(0.5*(pmax.x + pmin.x), pmax.y + 10 * (pmax.y - pmin.y),
			pmin.z - 100 * (pmax.z - pmin.z)));
		tetin.numberofpoints += 4;
		for (size_t i = N; i < N + 4; ++i)
		{
			tetin.pointlist[i * 3] = mesh_points_[i].x;
			tetin.pointlist[i * 3 + 1] = mesh_points_[i].y;
			tetin.pointlist[i * 3 + 2] = mesh_points_[i].z;
		}
	}
	// Run tetgen
	if (!voronoi)
		tetrahedralize("nQT1e-17", &tetin, &tetout);
	else
		tetrahedralize("nQvT1e-17", &tetin, &tetout);
}

vector<Vector3D> Voronoi3D::CreateBoundaryPoints(vector<vector<size_t> > const& to_duplicate,
	Vector3D const& ll, Vector3D const& ur)const
{
	vector<Face> faces = BuildBox(ll, ur);
	vector<Vector3D> res;
	for (size_t i = 0; i < to_duplicate.size(); ++i)
	{
		for (size_t j = 0; j < to_duplicate[i].size(); ++j)
			res.push_back(MirrorPoint(faces[i], mesh_points_[to_duplicate[i][j]]));
	}
	return res;
}

void Voronoi3D::Build(vector<Vector3D> const & points)
{
	Norg_ = points.size();
	// Clear data
	mesh_points_.clear();
	tetras_.clear();
	PointTetras_.clear();
	R_.clear();
	tetra_centers_.clear();
	// Voronoi Data
	FacesInCell_.clear();
	PointsInFace_.clear();
	FaceNeighbors_.clear();
	CM_.clear();
	volume_.clear();
	area_.clear();
	FacePoints_.clear();

	RunTetGen(points);
	CopyData();
	duplicated_points_ = FindBoundaryCandidates(ll_, ur_);
	vector<Vector3D> extra_points = CreateBoundaryPoints(duplicated_points_, ll_, ur_);
	// This could be made faster, also make sure no files are written
	mesh_points_.insert(mesh_points_.end(), extra_points.begin(), extra_points.end());
	RunTetGen(mesh_points_, true);

	CM_.resize(Norg_);
	volume_.resize(Norg_);
	// Copy the voronoi data
	CopyDataVoronoi();

	Norg_ = Norg_;
}

void Voronoi3D::CopyDataVoronoi()
{
	// copy tetra info
	size_t Ntetra = static_cast<size_t>(tetout.numberoftetrahedra);
	PointTetras_.clear();
	PointTetras_.resize(Norg_);
	for (size_t i = 0; i < Norg_; ++i)
		PointTetras_[i].reserve(20);
	for (size_t i = 0; i < Ntetra; ++i)
	{

		for (size_t j = 0; j < 4; ++j)
		{
			if (tetout.tetrahedronlist[i * 4 + j] < Norg_)
				PointTetras_[tetout.tetrahedronlist[i * 4 + j]].push_back(i);
		}
	}
	// calc all tetra radii
	R_.resize(Ntetra);
	boost::array<Vector3D, 4> points;
	for (size_t i = 0; i < Ntetra; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
			points[j] = mesh_points_[tetout.tetrahedronlist[i * 4 + j]];
		R_[i] = CalcRadius(points);
	}

	size_t Nface_points = static_cast<size_t>(tetout.numberofvpoints);
	FacePoints_.resize(Nface_points);
	for (size_t i = 0; i < Nface_points; ++i)
		FacePoints_[i] = Vector3D(tetout.vpointlist[3 * i], tetout.vpointlist[3 * i + 1], tetout.vpointlist[3 * i + 2]);
	FacesInCell_.clear();
	FacesInCell_.resize(Norg_);
	PointsInFace_.clear();
	FaceNeighbors_.clear();
	area_.clear();
	size_t Nfaces = static_cast<size_t>(tetout.numberofvfacets);
	area_.reserve(Nfaces);
	PointsInFace_.reserve(Nfaces);
	FaceNeighbors_.reserve(Nfaces);
	VecCompare compare(FacePoints_);
	vector<size_t> temp,temp2;
	for (size_t i = 0; i < Nfaces; ++i)
	{
		size_t N0 = tetout.vfacetlist[i].c1;
		size_t N1 = tetout.vfacetlist[i].c2;
		if (N0 >= Norg_&& N1 >= Norg_)
			continue;
		temp.clear();
		int nedges = tetout.vfacetlist[i].elist[0];
		for (int j = 0; j < nedges; ++j)
		{
			int edge = tetout.vfacetlist[i].elist[j + 1];
			temp.push_back(tetout.vedgelist[edge].v1);
			temp.push_back(tetout.vedgelist[edge].v2);
		}
		std::sort(temp.begin(), temp.end());
		temp = unique(temp);
		double R = 0;
		if (N0 < Norg_ && N1 < Norg_)
			R = std::min(R_[N0], R_[N1]);
		else
			R = N0 < Norg_ ? R_[N0] : R_[N1];
		RemoveDuplicates(FacePoints_, temp, R, temp2, compare);
		if (temp.size() >= 3)
		{
				ConvexHull3D(FacePoints_, temp, temp2);
				PointsInFace_.push_back(temp);
				FaceNeighbors_.push_back(std::pair<size_t, size_t>(N0, N1));
				area_.push_back(CalcFaceArea(temp, FacePoints_));
				if (N0 < Norg_)
					FacesInCell_[N0].push_back(area_.size() - 1);
				if (N1 < Norg_)
					FacesInCell_[N1].push_back(area_.size() - 1);
		}
	}
	for (size_t i = 0; i < Norg_; ++i)
		CalcCellCMVolume(i);
}

void Voronoi3D::CopyData()
{
	tetras_.resize(static_cast<size_t>(tetout.numberoftetrahedra));
	size_t Ntetra = tetras_.size();
	int counter = 0;
	PointTetras_.resize(Norg_ + 4);
	for (size_t i = 0; i < Ntetra; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			tetras_[i].points[j] = tetout.tetrahedronlist[counter];
			if (tetras_[i].points[j] >= Norg_)
				bigtet_ = i;
			tetras_[i].neighbors[j] = tetout.neighborlist[counter];
			++counter;
			PointTetras_[tetras_[i].points[j]].push_back(i);
		}
	}

	R_.resize(Ntetra, -1);
	tetra_centers_.resize(Ntetra);
}

void Voronoi3D::FillPointTetra(size_t point, size_t initetra)
{
	set_temp_.clear();
	for (size_t i = 0; i < 4; ++i)
	{
		int candidate = tetras_[initetra].neighbors[i];
		if (candidate >= 0)
			if (PointInVertices(tetras_[candidate].points, point) && set_temp_.find(candidate) == set_temp_.end())
			{
				set_temp_.insert(candidate);
				stack_temp_.push(candidate);
			}
	}
	while (!stack_temp_.empty())
	{
		int tocheck = stack_temp_.top();
		stack_temp_.pop();
		for (size_t i = 0; i < 4; ++i)
		{
			int candidate = tetras_[tocheck].neighbors[i];
			if (candidate >= 0)
				if (PointInVertices(tetras_[candidate].points, point) && set_temp_.find(candidate) == set_temp_.end())
				{
					set_temp_.insert(candidate);
					stack_temp_.push(candidate);
				}
		}
	}
	PointTetras_[point].resize(set_temp_.size());
	std::copy(set_temp_.begin(), set_temp_.end(), PointTetras_[point].begin());
}

double Voronoi3D::GetRadius(size_t index)
{
	if (R_[index] < 0)
		R_[index] = CalcTetraRadiusCenter(index);
	return R_[index];
}

double Voronoi3D::GetMaxRadius(size_t index)
{
	size_t N = PointTetras_[index].size();
	double res = 0;
	for (size_t i = 0; i < N; ++i)
		res = std::max(res, GetRadius(PointTetras_[index][i]));
	return res;
}

vector<vector<size_t> > Voronoi3D::FindIntersections(vector<Face>	const &box)
{
	size_t cur_loc = bigtet_;
	std::stack<size_t > check_stack;
	check_stack.push(cur_loc);
	size_t nbox = box.size();
	vector<vector<size_t> > res(nbox);
	Sphere sphere;
	vector<bool> checked(Norg_, false);
	while (!check_stack.empty())
	{
		cur_loc = check_stack.top();
		check_stack.pop();
		// Does tetra have any intersections?
		for (size_t k = 0; k < 4; ++k)
		{
			size_t point_tocheck = static_cast<size_t>(tetras_[cur_loc].points[k]);
			if (point_tocheck >= Norg_ || checked[point_tocheck])
				continue;
			checked[point_tocheck] = true;
			bool added = false;
			size_t Nneigh = PointTetras_[point_tocheck].size();

			for (size_t j = 0; j < nbox; ++j)
			{
				for (size_t l = 0; l < Nneigh; ++l)
				{
					size_t tetra_tocheck = PointTetras_[point_tocheck][l];
					sphere.radius = GetRadius(tetra_tocheck);
					sphere.center = tetra_centers_[tetra_tocheck];
					if (FaceSphereIntersections(box[j], sphere))
					{
						added = true;
						res[j].push_back(point_tocheck);
						break;
					}
				}
			}
			if (added)
				for (size_t j = 0; j < Nneigh; ++j)
					check_stack.push(PointTetras_[point_tocheck][j]);
		}
	}
	return res;
}

vector<vector<size_t> > Voronoi3D::FindBoundaryCandidates(Vector3D const& ll, Vector3D const& ur)
{
	// Build box polyhedra
	vector<Face> box = BuildBox(ll, ur);
	// run recursive script checking intersections, if small number use brute force
	if (Norg_ < 350)
	{
		vector<vector<size_t> > res(box.size());
		for (size_t i = 0; i < res.size(); ++i)
		{
			res[i].resize(Norg_);
			for (size_t j = 0; j < Norg_; ++j)
				res[i][j] = j;
		}
		return res;
	}
	vector<vector<size_t> > res = FindIntersections(box);
	// Clean up the results
	for (size_t i = 0; i < res.size(); ++i)
	{
		std::sort(res[i].begin(), res[i].end());
		res[i] = unique(res[i]);
	}
	return res;
}

double Voronoi3D::CalcTetraRadiusCenter(size_t index)
{
#define v1 (mesh_points_[tetras_[index].points[0]])
#define v2 (mesh_points_[tetras_[index].points[1]])
#define v3 (mesh_points_[tetras_[index].points[2]])
#define v4 (mesh_points_[tetras_[index].points[3]])

	Mat44<double> m_a(v1.x, v1.y, v1.z, 1,
		v2.x, v2.y, v2.z, 1,
		v3.x, v3.y, v3.z, 1,
		v4.x, v4.y, v4.z, 1);
	double a = m_a.determinant();

	Mat44<double> m_Dx(ScalarProd(v1, v1), v1.y, v1.z, 1,
		ScalarProd(v2, v2), v2.y, v2.z, 1,
		ScalarProd(v3, v3), v3.y, v3.z, 1,
		ScalarProd(v4, v4), v4.y, v4.z, 1);
	double Dx = m_Dx.determinant();

	Mat44<double> m_Dy(ScalarProd(v1, v1), v1.x, v1.z, 1,
		ScalarProd(v2, v2), v2.x, v2.z, 1,
		ScalarProd(v3, v3), v3.x, v3.z, 1,
		ScalarProd(v4, v4), v4.x, v4.z, 1);
	double Dy = -m_Dy.determinant();

	Mat44<double> m_Dz(ScalarProd(v1, v1), v1.x, v1.y, 1,
		ScalarProd(v2, v2), v2.x, v2.y, 1,
		ScalarProd(v3, v3), v3.x, v3.y, 1,
		ScalarProd(v4, v4), v4.x, v4.y, 1);
	double Dz = m_Dz.determinant();

	Mat44<double> m_c(ScalarProd(v1, v1), v1.x, v1.y, v1.z,
		ScalarProd(v2, v2), v2.x, v2.y, v2.z,
		ScalarProd(v3, v3), v3.x, v3.y, v3.z,
		ScalarProd(v4, v4), v4.x, v4.y, v4.z);
	double c = m_c.determinant();

#undef v1
#undef v2
#undef v3
#undef v4

	tetra_centers_[index] = Vector3D(Dx / (2 * a), Dy / (2 * a), Dz / (2 * a));
	return 0.5*sqrt(Dx*Dx + Dy*Dy + Dz*Dz - 4 * a*c) / std::abs(a);
}

Vector3D Voronoi3D::GetTetraCM(boost::array<Vector3D, 4> const& points)const
{
	Vector3D res;
	for (size_t i = 0; i < 4; ++i)
		res += 0.25*points[i];
	return res;
}

double Voronoi3D::GetTetraVolume(boost::array<Vector3D, 4> const& points)const
{
	Mat44<double> mat(points[0].x, points[0].y, points[0].z, 1,
		points[1].x, points[1].y, points[1].z, 1,
		points[2].x, points[2].y, points[2].z, 1,
		points[3].x, points[3].y, points[3].z, 1);
	double det = mat.determinant();
	return det / 6.0;
}

void Voronoi3D::CalcCellCMVolume(size_t index)
{
	// Make sure face is convexhull

	volume_[index] = 0;
	CM_[index] = Vector3D();
	size_t Nfaces = FacesInCell_[index].size();
	boost::array<Vector3D, 4> tetra;
	tetra[3] = mesh_points_[index];
	for (size_t i = 0; i < Nfaces; ++i)
	{
		size_t face = FacesInCell_[index][i];
		size_t Npoints = PointsInFace_[face].size();
		tetra[0] = FacePoints_[PointsInFace_[face][0]];

		for (size_t j = 0; j < Npoints - 2; ++j)
		{
			tetra[1] = FacePoints_[PointsInFace_[face][j + 1]];
			tetra[2] = FacePoints_[PointsInFace_[face][j + 2]];
			double vol = GetTetraVolume(tetra);
			volume_[index] += abs(vol);
			CM_[index] += vol*GetTetraCM(tetra);
		}
	}
	CM_[index] = CM_[index] / volume_[index];
}

namespace
{

	void binary_write_single_int(int n, std::ofstream& fh)
	{
		fh.write(reinterpret_cast<const char*>(&n), sizeof(int));
	}

	void binary_write_single_double(double d, std::ofstream& fh)
	{
		fh.write(reinterpret_cast<const char*>(&d), sizeof(double));
	}

	void binary_write_single_size_t(size_t n, std::ofstream& fh)
	{
		fh.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
	}
}

void Voronoi3D::output(std::string const& filename)const
{
	std::ofstream file_handle(filename.c_str(), std::ios::binary);

	//	binary_write_single_size_t(Norg_,file_handle);
	binary_write_single_int(Norg_, file_handle);
	//	binary_write_single_size_t(PointsInFace_.size(), file_handle); // Number of faces
	//	binary_write_single_size_t(FacePoints_.size(), file_handle); // Number of voronoi vertices points

		// Points
	for (size_t i = 0; i < Norg_; ++i)
	{
		binary_write_single_double(mesh_points_[i].x, file_handle);
		binary_write_single_double(mesh_points_[i].y, file_handle);
		binary_write_single_double(mesh_points_[i].z, file_handle);
	}

	// Find relevant edges
	vector<size_t> rel_faces;
	for (size_t i = 0; i < Norg_; ++i)
		for (size_t j = 0; j < FacesInCell_[i].size(); ++j)
			rel_faces.push_back(FacesInCell_[i][j]);
	sort(rel_faces.begin(), rel_faces.end());
	rel_faces = unique(rel_faces);
	vector<size_t> rel_edges;
	for (size_t i = 0; i < rel_faces.size(); ++i)
		for (size_t j = 0; j < tetout.vfacetlist[rel_faces[i]].elist[0]; ++j)
			rel_edges.push_back(tetout.vfacetlist[rel_faces[i]].elist[j + 1]);
	sort(rel_edges.begin(), rel_edges.end());
	rel_edges = unique(rel_edges);


	int nedges = rel_edges.size();
	binary_write_single_int(nedges, file_handle);
	for (size_t i = 0; i < nedges; ++i)
	{
		binary_write_single_int(tetout.vedgelist[rel_edges[i]].v1, file_handle);
		binary_write_single_int(tetout.vedgelist[rel_edges[i]].v2, file_handle);
	}
	binary_write_single_int(FacePoints_.size(), file_handle);
	// Face Points
	for (size_t i = 0; i < FacePoints_.size(); ++i)
	{
		binary_write_single_double(FacePoints_[i].x, file_handle);
		binary_write_single_double(FacePoints_[i].y, file_handle);
		binary_write_single_double(FacePoints_[i].z, file_handle);
	}
	file_handle.close();
}

size_t Voronoi3D::GetPointNo(void) const
{
	return Norg_;
}

Vector3D Voronoi3D::GetMeshPoint(size_t index) const
{
	return mesh_points_[index];
}

double Voronoi3D::GetArea(size_t index) const
{
	return area_[index];
}

Vector3D const& Voronoi3D::GetCellCM(size_t index) const
{
	return CM_[index];
}

size_t Voronoi3D::GetTotalFacesNumber(void) const
{
	return FaceNeighbors_.size();
}

double Voronoi3D::GetWidth(size_t index) const
{
	return pow(3 * volume_[index] * 0.25 / M_PI, 0.333333333);
}

double Voronoi3D::GetVolume(size_t index) const
{
	return volume_[index];
}

vector<size_t>const& Voronoi3D::GetCellFaces(size_t index) const
{
	return FacesInCell_[index];
}

vector<Vector3D>& Voronoi3D::GetMeshPoints(void)
{
	return mesh_points_;
}

vector<size_t> Voronoi3D::GetNeighbors(size_t index)const
{
	size_t N = FacesInCell_[index].size();
	vector<size_t> res(N);
	for (size_t i = 0; i < N; ++i)
	{
		size_t face = FacesInCell_[index][i];
		res[i] = FaceNeighbors_[face].first == index ? FaceNeighbors_[face].second :
			FaceNeighbors_[face].first;
	}
	return res;
}

Tessellation3D* Voronoi3D::clone(void) const
{
	return new Voronoi3D(*this);
}

bool Voronoi3D::NearBoundary(size_t index) const
{
	size_t N = FacesInCell_[index].size();
	for (size_t i = 0; i < N; ++i)
	{
		if (BoundaryFace(FacesInCell_[index][i]))
			return true;
	}
	return false;
}

bool Voronoi3D::BoundaryFace(size_t index) const
{
	if (FaceNeighbors_[index].first >= Norg_ || FaceNeighbors_[index].second >= Norg_)
		return true;
	else
		return false;
}

vector<vector<size_t> >& Voronoi3D::GetDuplicatedPoints(void)
{
	return duplicated_points_;
}

vector<vector<size_t> >const& Voronoi3D::GetDuplicatedPoints(void)const
{
	return duplicated_points_;
}

size_t Voronoi3D::GetTotalPointNumber(void)const
{
	return mesh_points_.size();
}

vector<Vector3D>& Voronoi3D::GetAllCM(void)
{
	return CM_;
}

void Voronoi3D::GetNeighborNeighbors(vector<size_t> &result, size_t point)const
{
	result.clear();
	result.reserve(70);
	vector<size_t> neigh = GetNeighbors(point);
	result = neigh;
	size_t N = neigh.size();
	vector<size_t> temp;
	for (size_t i = 0; i < N; ++i)
	{
		if (neigh[i] < Norg_)
		{
			temp = GetNeighbors(neigh[i]);
			result.insert(result.end(), temp.begin(), temp.end());
		}
	}
	sort(result.begin(), result.end());
	result = unique(result);
}

Vector3D Voronoi3D::Normal(size_t faceindex)const
{
	return mesh_points_[FaceNeighbors_[faceindex].second] - mesh_points_[FaceNeighbors_[faceindex].first];
}

bool Voronoi3D::IsGhostPoint(size_t index)const
{
	return index >= Norg_;
}

Vector3D Voronoi3D::FaceCM(size_t index)const
{
	size_t N = PointsInFace_[index].size();
	Vector3D res = FacePoints_[PointsInFace_[index][0]];
	for (size_t i = 1; i < N; ++i)
		res += FacePoints_[PointsInFace_[index][i]];
	return res / N;
}

Vector3D Voronoi3D::CalcFaceVelocity(size_t index, Vector3D const& v0, Vector3D const& v1)const
{
	size_t p0 = FaceNeighbors_[index].first;
	size_t p1 = FaceNeighbors_[index].second;
	Vector3D r0 = GetMeshPoint(p0);
	Vector3D r1 = GetMeshPoint(p1);
	Vector3D r_diff = r1 - r0;
	double abs_r_diff = abs(r_diff);
	
	Vector3D f = FaceCM(index);

	Vector3D delta_w = ScalarProd((v0 - v1), (f - (r1+r0) / 2)) * r_diff / (abs_r_diff * abs_r_diff);
	Vector3D w = (v0 + v1) / 2 + delta_w;
	return w;
}

vector<Vector3D>const& Voronoi3D::GetFacePoints(void) const
{
	return FacePoints_;
}

vector<size_t>const& Voronoi3D::GetPointsInFace(size_t index) const
{
	return PointsInFace_[index];
}