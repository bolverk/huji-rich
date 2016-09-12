#include "Voronoi3D.hpp"
#include <algorithm>
#include <stack>
#include "Mat44.hpp"
#include "../../misc/utils.hpp"
#include <fstream>
#include <iostream>
#include <boost/container/flat_map.hpp>

bool PointInPoly(Tessellation3D const& tess, Vector3D const& point, size_t index)
{
	vector<size_t> const& faces = tess.GetCellFaces(index);
	size_t N = faces.size();
	for (size_t i = 0; i < N; ++i)
		if (ScalarProd(tess.GetMeshPoint(index) - tess.GetFacePoints()[tess.GetPointsInFace(faces[i])[0]],
			point - tess.GetFacePoints()[tess.GetPointsInFace(faces[i])[0]]) < 0)
			return false;
	return true;
}


namespace
{
	void BuildBigTetra(vector<Vector3D> & mesh_points)
	{
		Vector3D ll(mesh_points[0]), ur(mesh_points[0]);
		size_t N = mesh_points.size();
		for (size_t i = 0; i < N; ++i)
		{
			ll.x = std::min(ll.x, mesh_points[i].x);
			ll.y = std::min(ll.y, mesh_points[i].y);
			ll.z = std::min(ll.z, mesh_points[i].z);
			ur.x = std::max(ur.x, mesh_points[i].x);
			ur.y = std::max(ur.y, mesh_points[i].y);
			ur.z = std::max(ur.z, mesh_points[i].z);
		}
		double bigratio = 3;
		mesh_points.push_back(Vector3D(0.5*(ll.x + ur.x), 0.5*(ur.y + ll.y), ur.z + bigratio * (ur.z - ll.z)));
		mesh_points.push_back(Vector3D(ll.x - bigratio * (ur.x - ll.x), ll.y - bigratio * (ur.y - ll.y),
			ll.z - bigratio * (ur.z - ll.z)));
		mesh_points.push_back(Vector3D(ur.x + bigratio * (ur.x - ll.x), ll.y - bigratio * (ur.y - ll.y),
			ll.z - bigratio * (ur.z - ll.z)));
		mesh_points.push_back(Vector3D(0.5*(ur.x + ll.x), ur.y + bigratio * (ur.y - ll.y),
			ll.z - bigratio * (ur.z - ll.z)));
	}

#ifdef RICH_MPI
	void TalkSymmetry(vector<int> & to_talk_with)
	{
		int wsize;
		MPI_Comm_size(MPI_COMM_WORLD, &wsize);
		vector<int> totalk(static_cast<size_t>(wsize), 0);
		vector<int> scounts(totalk.size(), 1);
		for (size_t i = 0; i < to_talk_with.size(); ++i)
			totalk[to_talk_with[i]] = 1;
		int nrecv=0;
		MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
			MPI_COMM_WORLD);

		vector<MPI_Request> req(to_talk_with.size());
		for (size_t i = 0; i < to_talk_with.size(); ++i)
			MPI_Isend(&wsize, 1, MPI_INT, to_talk_with[i], 3, MPI_COMM_WORLD, &req[i]);
		vector<int> talkwithme;
		for (int i = 0; i < nrecv; ++i)
		{
			MPI_Status status;
			MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
			talkwithme.push_back(status.MPI_SOURCE);
		}
		MPI_Waitall(static_cast<int>(to_talk_with.size()), &req[0], MPI_STATUSES_IGNORE);
		vector<int> new_talk_with_me;
		for (size_t i = 0; i < to_talk_with.size(); ++i)
			if (std::find(talkwithme.begin(), talkwithme.end(), to_talk_with[i]) != talkwithme.end())
				new_talk_with_me.push_back(to_talk_with[i]);
		to_talk_with = new_talk_with_me;
	}
#endif //RICH_MPI

	double CalcFaceArea(vector<size_t> const& indeces, vector<Vector3D> const& points)
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
			if (std::abs(points[i].x - points[j].x) > 1e-8*R)
				return points[i].x < points[j].x;
			else
				if (std::abs(points[i].y - points[j].y) > 1e-8*R)
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

	void ConvexHull3D(vector<Vector3D> const& points, vector<size_t> &indeces, vector<size_t> & temp,
		Vector3D const& normal)
	{
		// Find center point in plane
		Vector3D center = (points[indeces[0]] + points[indeces[1]] + points[indeces[2]]);
		for (size_t i = 3; i < indeces.size(); ++i)
			center = center + points[indeces[i]];
		center = center / static_cast<double>(indeces.size());
		vector<double> angles(indeces.size() - 1);
		size_t Nloop = angles.size();
		Vector3D main_vector = points[indeces[Nloop]] - center;
		main_vector = main_vector / abs(main_vector);
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

	bool PointInDomain(Vector3D const& ll, Vector3D const& ur, Vector3D const& point)
	{
		if (point.x > ll.x&&point.x<ur.x&&point.y>ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			return true;
		else
			return false;
	}

	Vector3D MirrorPoint(Face const& face, Vector3D const& point)
	{
		Vector3D normal = CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
		normal = normal / abs(normal);
		return point - (2 * ScalarProd(point - face.vertices[0], normal))*normal;
	}
}

#ifdef RICH_MPI
vector<Vector3D> Voronoi3D::UpdateMPIPoints(Tessellation3D const& vproc, int rank,
	vector<Vector3D> const& points,vector<size_t> &selfindex,vector<int> &sentproc,
	vector<vector<size_t> > &sentpoints)
{
	vector<Vector3D> res;
	res.reserve(points.size());
	selfindex.clear();
	size_t npoints = points.size();
	size_t nproc = vproc.GetPointNo();
	vector<size_t> neighbors = vproc.GetNeighbors(static_cast<size_t>(rank));
	vector<size_t> realneigh;
	sentpoints.clear();
	sentproc.clear();
	for (size_t i = 0; i < neighbors.size(); ++i)
		if (static_cast<size_t>(neighbors[i]) < nproc)
		{
			realneigh.push_back(neighbors[i]);
			sentproc.push_back(static_cast<int>(neighbors[i]));
		}
	size_t Nreal = realneigh.size();
	sentpoints.resize(sentproc.size());

	for (size_t i = 0; i<npoints; ++i)
	{
		Vector3D temp = points[i];
		if(PointInPoly(vproc,temp,static_cast<size_t>(rank)))
		{
			res.push_back(temp);
			selfindex.push_back(i);
			continue;
		}
		bool good = false;
		for (size_t j = 0; j<Nreal; ++j)
		{
			if (PointInPoly(vproc, temp, realneigh[j]))
			{
				sentpoints[j].push_back(i);
				good = true;
				break;
			}
		}
		if (good)
			continue;
		for (size_t j = 0; j<nproc; ++j)
		{
			if (std::find(realneigh.begin(), realneigh.end(), j) != realneigh.end() || j == static_cast<size_t>(rank))
				continue;
			if (PointInPoly(vproc, temp, j))
			{
				good = true;
				size_t index = std::find(sentproc.begin(), sentproc.end(), j) - sentproc.begin();
				if (index >= sentproc.size())
				{
					sentproc.push_back(static_cast<int>(j));
					sentpoints.push_back(vector<size_t>(1, i));
				}
				else
					sentpoints[index].push_back(i);
				break;
			}
		}
		if (good)
			continue;
		UniversalError eo("Point is not inside any processor");
		eo.AddEntry("CPU rank", rank);
		eo.AddEntry("Point number", static_cast<double>(i));
		eo.AddEntry("Point x cor", points[i].x);
		eo.AddEntry("Point y cor", points[i].y);
		throw eo;
	}
	// Send/Recv the points
	// Communication
	int wsize;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	vector<int> totalk(static_cast<size_t>(wsize), 0);
	vector<int> scounts(totalk.size(), 1);
	for (size_t i = 0; i < sentproc.size(); ++i)
		totalk[sentproc[i]] = 1;
	int nrecv;
	MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

	vector<MPI_Request> req(sentproc.size());
	for (size_t i = 0; i < sentproc.size(); ++i)
		MPI_Isend(&wsize, 1, MPI_INT, sentproc[i], 3, MPI_COMM_WORLD, &req[i]);
	vector<int> talkwithme;
	for (int i = 0; i < nrecv; ++i)
	{
		MPI_Status status;
		MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
		talkwithme.push_back(status.MPI_SOURCE);
	}
	MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	for (size_t i = 0; i < talkwithme.size(); ++i)
	{
		if (std::find(sentproc.begin(), sentproc.end(), talkwithme[i]) == sentproc.end())
		{
			sentproc.push_back(talkwithme[i]);
			sentpoints.push_back(vector<size_t>());
		}
	}
	// Point exchange
	vector<vector<Vector3D> > incoming=MPI_exchange_data(sentproc, sentpoints, points);
	// Combine the vectors
	for (size_t i = 0; i < incoming.size(); ++i)
		for (size_t j = 0; j < incoming[i].size(); ++j)
			res.push_back(incoming[i][j]);
	return res;
}
#endif //RICH_MPI



Voronoi3D::Voronoi3D()
{}

Voronoi3D::Voronoi3D(Vector3D const& ll, Vector3D const& ur) :ll_(ll), ur_(ur) {}

Tetrahedron::Tetrahedron(void) {}

void Voronoi3D::RunTetGen(vector<Vector3D> const& points, tetgenio &tetin, tetgenio &tetout, bool voronoi)
{
	tetin.firstnumber = 0;
	tetin.numberofpoints = static_cast<int>(points.size());
	size_t N = points.size();
	tetin.pointlist = new double[N * 3];
	for (size_t i = 0; i < N; ++i)
	{
		tetin.pointlist[i * 3] = points[i].x;
		tetin.pointlist[i * 3 + 1] = points[i].y;
		tetin.pointlist[i * 3 + 2] = points[i].z;
	}
	// Run tetgen
	if (!voronoi)
	{
		char msg[] = "nQT1e-17";
		tetrahedralize(msg, &tetin, &tetout);
	}
	else
	{
		char msg[] = "nQvT1e-17";
		tetrahedralize(msg, &tetin, &tetout);
	}
}

void Voronoi3D::CalcRigidCM(size_t face_index)
{
	Vector3D normal = normalize(mesh_points_[FaceNeighbors_[face_index].first] - mesh_points_[FaceNeighbors_[face_index].second]);
	size_t real, other;
	if (FaceNeighbors_[face_index].first >= Norg_)
	{
		real = FaceNeighbors_[face_index].second;
		other = FaceNeighbors_[face_index].first;
	}
	else
	{
		real = FaceNeighbors_[face_index].first;
		other = FaceNeighbors_[face_index].second;
	}
	CM_[other] = CM_[real] - 2 * normal*ScalarProd(normal, CM_[real]-FacePoints_[PointsInFace_[face_index][0]]);
}



vector<Vector3D> Voronoi3D::CreateBoundaryPoints(vector<std::pair<size_t,size_t> > const& to_duplicate)
{
	self_duplicate_.resize(to_duplicate.size());
	vector<Face> faces = BuildBox(ll_, ur_);
	vector<Vector3D> res;
	for (size_t i = 0; i < to_duplicate.size(); ++i)
	{
		res.push_back(MirrorPoint(faces[to_duplicate[i].first], mesh_points_[to_duplicate[i].second]));
		self_duplicate_[i] = to_duplicate[i].second;
	}
	return res;
}

#ifdef RICH_MPI
vector<Vector3D> Voronoi3D::CreateBoundaryPointsMPI(vector<std::pair<size_t, size_t> > const& to_duplicate,
	Tessellation3D const& tproc)
{
	boost::container::flat_map<size_t, Face> boundary_faces;
	self_duplicate_.clear();
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::set<int> totalkwith;
	size_t Nproc = tproc.GetPointNo();
	for (size_t i = 0; i < to_duplicate.size(); ++i)
	{
		std::pair<size_t, size_t>const& neigh = tproc.GetFaceNeighbors(to_duplicate[i].first);
		size_t other = (neigh.first == static_cast<size_t> (rank)) ? neigh.second : neigh.first;
		if (other < Nproc)
			totalkwith.insert(static_cast<int>(other));
	}
	duplicatedprocs_.resize(totalkwith.size());
	size_t counter = 0;
	for (std::set<int>::iterator it = totalkwith.begin(); it != totalkwith.end(); ++it,++counter)
	{
		duplicatedprocs_[counter] = *it;
	}
	TalkSymmetry(duplicatedprocs_);
	duplicated_points_.resize(duplicatedprocs_.size());
	std::fill(duplicated_points_.begin(), duplicated_points_.end(), vector<size_t>());
	vector<Vector3D> res;
	// Get the indeces and deal with selfboundary
	for (size_t i = 0; i < to_duplicate.size(); ++i)
	{
		std::pair<size_t, size_t>const& neigh = tproc.GetFaceNeighbors(to_duplicate[i].first);
		size_t other = (neigh.first == static_cast<size_t> (rank)) ? neigh.second : neigh.first;
		if (other < Nproc)
		{
			size_t index = std::lower_bound(duplicatedprocs_.begin(), duplicatedprocs_.end(), static_cast<int>(other))
				- duplicatedprocs_.begin();
			if(index<duplicatedprocs_.size())
				duplicated_points_[index].push_back(neigh.second);
		}
		else
		{
			if (boundary_faces.find(neigh.first) == boundary_faces.end())
			{
				Face f(VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(neigh.first)), 0, 0);
				boundary_faces.insert(std::pair<size_t, Face>(neigh.first, f));
			}
			self_duplicate_.push_back(neigh.second);
			res.push_back(MirrorPoint(boundary_faces[neigh.first], mesh_points_[neigh.second]));
		}
	}
	// Communicate
	vector<vector<Vector3D> > toadd = MPI_exchange_data(duplicatedprocs_, duplicated_points_, mesh_points_);
	// Add points
	Nghost_.clear();
	Nghost_.resize(toadd.size());
	for (size_t i = 0; i < toadd.size(); ++i)
		for (size_t j = 0; j < toadd[i].size(); ++j)
		{
			Nghost_[i].push_back(Norg_ + 4 + res.size());
			res.push_back(toadd[i][j]);
		}
	return res;
}
#endif //RICH_MPI


vector<vector<size_t> > const& Voronoi3D::GetGhostIndeces(void) const
{
	return Nghost_;
}

#ifdef RICH_MPI
void Voronoi3D::Build(vector<Vector3D> const & points,Tessellation3D const& tproc)
{
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

	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Vector3D> new_points = UpdateMPIPoints(tproc, rank,points,self_index_,sentprocs_,sentpoints_);
	mesh_points_.reserve(new_points.size() + 4);
	mesh_points_ = new_points;
	Norg_ = mesh_points_.size();
	// Build big tetrahedron
	BuildBigTetra(mesh_points_);

	tetgenio tetin, tetout;
	RunTetGen(mesh_points_,tetin,tetout);
	CopyData(tetout);

	vector<std::pair<size_t, size_t> > ghost_index = FindIntersections(tproc,false);
	vector<Vector3D> extra_points = CreateBoundaryPointsMPI(ghost_index,tproc);
	mesh_points_.insert(mesh_points_.end(), extra_points.begin(), extra_points.end());
	tetin.deinitialize();
	tetin.initialize();
	tetout.deinitialize();
	tetout.initialize();
	RunTetGen(mesh_points_, tetin, tetout);
	ghost_index = FindIntersections(tproc, true);
	extra_points = CreateBoundaryPointsMPI(ghost_index, tproc);
	mesh_points_.resize(Norg_ + 4);

	// This could be made faster, also make sure no files are written
	mesh_points_.insert(mesh_points_.end(), extra_points.begin(), extra_points.end());

	tetin.deinitialize();
	tetin.initialize();
	tetout.deinitialize(); 
	tetout.initialize();
	RunTetGen(mesh_points_,tetin,tetout,true);

	CM_.resize(mesh_points_.size());
	volume_.resize(Norg_);
	// Copy the voronoi data
	CopyDataVoronoi(tetout);

	Norg_ = Norg_;
}
#endif


void Voronoi3D::Build(vector<Vector3D> const & points)
{
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
	mesh_points_.reserve(points.size() + 4);
	mesh_points_ = points;
	Norg_ = mesh_points_.size();
	// Build big tetrahedron
	BuildBigTetra(mesh_points_);

	tetgenio tetin, tetout;
	RunTetGen(mesh_points_, tetin, tetout);
	CopyData(tetout);

	vector<std::pair<size_t, size_t> > ghost_index = SerialFindIntersections();
	vector<Vector3D> extra_points = CreateBoundaryPoints(ghost_index);

	// This could be made faster, also make sure no files are written
	mesh_points_.insert(mesh_points_.end(), extra_points.begin(), extra_points.end());

	tetin.deinitialize();
	tetin.initialize();
	tetout.deinitialize();
	tetout.initialize();
	RunTetGen(mesh_points_, tetin, tetout, true);

	CM_.resize(mesh_points_.size());
	volume_.resize(Norg_);
	// Copy the voronoi data
	CopyDataVoronoi(tetout);

	Norg_ = Norg_;
}



void Voronoi3D::CopyDataVoronoi(tetgenio &tetin)
{
	// copy tetra info
	size_t Ntetra = static_cast<size_t>(tetin.numberoftetrahedra);
	PointTetras_.clear();
	PointTetras_.resize(Norg_);
	for (size_t i = 0; i < Norg_; ++i)
		PointTetras_[i].reserve(30);
	for (size_t i = 0; i < Ntetra; ++i)
	{

		for (size_t j = 0; j < 4; ++j)
		{
			if (tetin.tetrahedronlist[i * 4 + j] < static_cast<int>(Norg_))
				PointTetras_[tetin.tetrahedronlist[i * 4 + j]].push_back(i);
		}
	}
	// calc all tetra radii
	R_.resize(Ntetra);
	boost::array<Vector3D, 4> points;
	for (size_t i = 0; i < Ntetra; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
			points[j] = mesh_points_[tetin.tetrahedronlist[i * 4 + j]];
		R_[i] = CalcRadius(points);
	}

	size_t Nface_points = static_cast<size_t>(tetin.numberofvpoints);
	FacePoints_.resize(Nface_points);
	for (size_t i = 0; i < Nface_points; ++i)
		FacePoints_[i] = Vector3D(tetin.vpointlist[3 * i], tetin.vpointlist[3 * i + 1], tetin.vpointlist[3 * i + 2]);
	FacesInCell_.clear();
	FacesInCell_.resize(Norg_);
	PointsInFace_.clear();
	FaceNeighbors_.clear();
	area_.clear();
	size_t Nfaces = static_cast<size_t>(tetin.numberofvfacets);
	area_.reserve(Nfaces);
	PointsInFace_.reserve(Nfaces);
	FaceNeighbors_.reserve(Nfaces);
	VecCompare compare(FacePoints_);
	vector<size_t> temp, temp2;
	vector<size_t> boundary_faces;
	for (size_t i = 0; i < Nfaces; ++i)
	{
		size_t N0 = tetin.vfacetlist[i].c1;
		size_t N1 = tetin.vfacetlist[i].c2;
		if (N0 >= Norg_&& N1 >= Norg_)
			continue;
		temp.clear();
		int nedges = tetin.vfacetlist[i].elist[0];
		for (int j = 0; j < nedges; ++j)
		{
			int edge = tetin.vfacetlist[i].elist[j + 1];
			if (edge >= 0)
			{
				int etemp = tetin.vedgelist[edge].v1;
				if (etemp >= 0)
					temp.push_back(tetin.vedgelist[edge].v1);
				etemp = tetin.vedgelist[edge].v2;
				if (etemp >= 0)
					temp.push_back(tetin.vedgelist[edge].v2);
			}
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
			Vector3D normal = normalize(mesh_points_[N0] - mesh_points_[N1]);
			ConvexHull3D(FacePoints_, temp, temp2, normal);
			PointsInFace_.push_back(temp);
			FaceNeighbors_.push_back(std::pair<size_t, size_t>(N0, N1));
			area_.push_back(CalcFaceArea(temp, FacePoints_));
			if (N0 < Norg_)
				FacesInCell_[N0].push_back(area_.size() - 1);
			if (N1 < Norg_)
				FacesInCell_[N1].push_back(area_.size() - 1);
			if (N1 >= Norg_ || N0 >= Norg_)
				boundary_faces.push_back(area_.size() - 1);
		}
	}
	for (size_t i = 0; i < Norg_; ++i)
		CalcCellCMVolume(i);
	size_t N = boundary_faces.size();
	for (size_t i = 0; i < N; ++i)
		CalcRigidCM(i);
}

void Voronoi3D::CopyData(tetgenio &tetout)
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
			if (tetras_[i].points[j] >= static_cast<int>(Norg_))
				bigtet_ = i;
			tetras_[i].neighbors[j] = tetout.neighborlist[counter];
			++counter;
			PointTetras_[tetras_[i].points[j]].push_back(i);
		}
	}

	R_.resize(Ntetra);
	std::fill(R_.begin(), R_.end(), -1);
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
	return 2*res;
}

vector<size_t>  Voronoi3D::FindIntersectionsSingle(vector<Face> const& box, size_t point, Sphere &sphere)
{
	size_t N = PointTetras_[point].size();
	vector<size_t> res;
	res.reserve(4);
	for (size_t j = 0; j < box.size(); ++j)
	{
		for (size_t i = 0; i < N; ++i)
		{
			sphere.radius = GetRadius(PointTetras_[point][i]);
			sphere.center = tetra_centers_[PointTetras_[point][i]];
			if (FaceSphereIntersections(box[j], sphere))
			{
				res.push_back(j);
				break;
			}
		}
	}
	return res;
}

vector<size_t>  Voronoi3D::FindIntersectionsRecursive(Tessellation3D const& tproc,size_t rank,size_t point,
	Sphere &sphere,bool recursive)
{
	vector<size_t> res;
	size_t N = tproc.GetPointNo();
	size_t Nfaces = tproc.GetTotalFacesNumber();
	vector<bool> visited(Nfaces, false);
	std::stack<size_t> to_check;
	size_t Ntetra=PointTetras_[point].size();
	for (size_t j = 0; j < Ntetra; ++j)
	{
		sphere.radius = GetRadius(PointTetras_[point][j]);
		sphere.center = tetra_centers_[PointTetras_[point][j]];
		vector<size_t> faces = tproc.GetCellFaces(rank);
		for (size_t i = 0; i < faces.size(); ++i)
			to_check.push(faces[i]);
		while (!to_check.empty())
		{
			size_t cur = to_check.top();
			to_check.pop();
			if (visited[cur])
				continue;
			visited[cur] = true;
			Face f(VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(cur)), tproc.GetFaceNeighbors(cur).first,
				tproc.GetFaceNeighbors(cur).second);
			if (FaceSphereIntersections(f, sphere))
			{
				res.push_back(cur);
				if (recursive)
				{
					if (f.neighbors.first < N)
					{
						vector<size_t> const& faces_temp = tproc.GetCellFaces(f.neighbors.first);
						for (size_t i = 0; i < faces_temp.size(); ++i)
							if (!visited[faces_temp[i]])
								to_check.push(faces_temp[i]);
					}
					if (f.neighbors.second < N)
					{
						vector<size_t> const& faces_temp = tproc.GetCellFaces(f.neighbors.second);
						for (size_t i = 0; i < faces_temp.size(); ++i)
							if (!visited[faces_temp[i]])
								to_check.push(faces_temp[i]);
					}
				}
			}
		}
	}
	std::sort(res.begin(), res.end());
	return unique(res);
}


void Voronoi3D::GetPointToCheck(size_t point, vector<bool> const& checked,vector<size_t> &res)
{
	res.clear();
	size_t ntetra=PointTetras_[point].size();
	for (size_t i = 0; i < ntetra; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
			if (tetras_[PointTetras_[point][i]].points[j]<static_cast<int>(Norg_)&&!checked[tetras_[PointTetras_[point][i]].points[j]])
				res.push_back(tetras_[PointTetras_[point][i]].points[j]);
	}
	std::sort(res.begin(), res.end());
	res = unique(res);
}

size_t Voronoi3D::GetFirstPointToCheck(void)const
{
	for (size_t i = 0; i < 4; ++i)
		if (tetras_[bigtet_].points[i] < static_cast<int>(Norg_))
			return tetras_[bigtet_].points[i];
	throw UniversalError("Can't find first point to start boundary search");
}

#ifdef RICH_MPI
vector<std::pair<size_t,size_t> > Voronoi3D::FindIntersections(Tessellation3D const& tproc,bool recursive)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<Face> box = BuildBox(ll_, ur_);
	size_t cur_loc = GetFirstPointToCheck();
	std::stack<size_t > check_stack;
	vector<size_t> point_neigh;
	check_stack.push(cur_loc);
	//size_t nbox = box.size();
	vector<std::pair<size_t, size_t> > res;
	Sphere sphere;
	vector<bool> checked(Norg_, false),will_check(Norg_,false);
	will_check[cur_loc] = true;
	while (!check_stack.empty())
	{
		cur_loc = check_stack.top();
		check_stack.pop();
		checked[cur_loc] = true;
		// Does sphere have any intersections?
		bool added = false;
		vector<size_t> intersecting_faces = FindIntersectionsRecursive(tproc, static_cast<size_t>(rank),cur_loc, sphere,recursive);
		if (!intersecting_faces.empty())
		{
			added = true;
			for (size_t j = 0; j < intersecting_faces.size(); ++j)
				res.push_back(std::pair<size_t, size_t>(intersecting_faces[j], cur_loc));
		}
		if (added)
		{
			GetPointToCheck(cur_loc, checked, point_neigh);
			size_t Nneigh = point_neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
				if (point_neigh[j] < Norg_ && !will_check[point_neigh[j]])
				{
					check_stack.push(point_neigh[j]);
					will_check[point_neigh[j]] = true;
				}
		}
	}
	return res;
}
#endif //RICH_MPI

vector<std::pair<size_t, size_t> > Voronoi3D::SerialFindIntersections()
{
	vector<Face> box = BuildBox(ll_, ur_);
	size_t cur_loc = GetFirstPointToCheck();
	std::stack<size_t > check_stack;
	vector<size_t> point_neigh;
	check_stack.push(cur_loc);
	//size_t nbox = box.size();
	vector<std::pair<size_t, size_t> > res;
	Sphere sphere;
	vector<bool> checked(Norg_, false), will_check(Norg_, false);
	will_check[cur_loc] = true;
	while (!check_stack.empty())
	{
		cur_loc = check_stack.top();
		check_stack.pop();
		checked[cur_loc] = true;
		// Does sphere have any intersections?
		bool added = false;
		vector<size_t> intersecting_faces = FindIntersectionsSingle(box, cur_loc, sphere);
		if (!intersecting_faces.empty())
		{
			added = true;
			for (size_t j = 0; j < intersecting_faces.size(); ++j)
				res.push_back(std::pair<size_t, size_t>(intersecting_faces[j], cur_loc));
		}
		if (added)
		{
			GetPointToCheck(cur_loc, checked, point_neigh);
			size_t Nneigh = point_neigh.size();
			for (size_t j = 0; j < Nneigh; ++j)
				if (point_neigh[j] < Norg_ && !will_check[point_neigh[j]])
				{
					check_stack.push(point_neigh[j]);
					will_check[point_neigh[j]] = true;
				}
		}
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
			volume_[index] += std::abs(vol);
			CM_[index] += std::abs(vol)*GetTetraCM(tetra);
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

void Voronoi3D::output(std::string const& /*filename*/)const
{
	/*
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
	file_handle.close();*/
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
	{
#ifdef RICH_MPI
		if (PointInDomain(ll_, ur_, mesh_points_[index]))
			return false;
		else
#endif
			return true;
	}
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
	return res / static_cast<double>(N);
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

	Vector3D delta_w = ScalarProd((v0 - v1), (f - (r1 + r0) / 2)) * r_diff / (abs_r_diff * abs_r_diff);
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

std::pair<size_t, size_t> Voronoi3D::GetFaceNeighbors(size_t face_index)const
{
	return std::pair<size_t, size_t>(FaceNeighbors_[face_index]);
}

vector<int> Voronoi3D::GetDuplicatedProcs(void)const
{
	return duplicatedprocs_;
}

vector<int> Voronoi3D::GetSentProcs(void)const
{
	return sentprocs_;
}

vector<vector<size_t> > const& Voronoi3D::GetSentPoints(void)const
{
	return sentpoints_;
}

vector<size_t> const& Voronoi3D::GetSelfIndex(void) const
{
	return self_index_;
}

vector<size_t> const& Voronoi3D::GetSelfDuplicate(void)const
{
	return self_duplicate_;
}
