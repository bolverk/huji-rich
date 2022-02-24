#include "mesh_generator3D.hpp"
#include <array>
#ifdef RICH_MPI
#include <mpi.h>
#endif

using std::array;

 vector<Vector3D> CartesianMesh(std::size_t nx, std::size_t ny, std::size_t nz, Vector3D const& lower_left, Vector3D const& upper_right,
	Voronoi3D const* tproc)
{
	assert(upper_right.x > lower_left.x);
	assert(upper_right.y > lower_left.y);
	assert(upper_right.z > lower_left.z);

	vector<Vector3D> res;
	const double dx = (upper_right.x - lower_left.x) /
		static_cast<double>(nx);
	const double dy = (upper_right.y - lower_left.y) /
		static_cast<double>(ny);
	const double dz = (upper_right.z - lower_left.z) /
		static_cast<double>(nz);
	int rank = 0;
#ifdef RICH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	for(size_t i = 0; i < nx; ++i)
		for(size_t j = 0; j < ny; ++j)
			for(size_t k = 0; k < nz; ++k)
			{

				Vector3D new_point = Vector3D(lower_left.x+0.5*dx+i*dx, lower_left.y+0.5*dy+j*dy, lower_left.z+0.5*dz+k*dz);
				if (tproc == nullptr)
					res.push_back(new_point);
				else
					if (PointInPoly(*tproc, new_point, static_cast<size_t>(rank)))
						res.push_back(new_point);	
			}
	return res;
}

vector<Vector3D> RandRectangular(std::size_t PointNum, Vector3D const& ll, 
	Vector3D const& ur, Voronoi3D const* tproc)
{
	typedef boost::mt19937_64 base_generator_type;
	double ran[3];
	Vector3D diff = ur - ll;
	vector<Vector3D> res;
	Vector3D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	if (tproc == nullptr)
	{
		for (size_t i = 0; i < PointNum; ++i)
		{
			ran[0] = dist(generator);
			ran[1] = dist(generator);
			ran[2] = dist(generator);
			point.x = ran[0] * diff.x + ll.x;
			point.y = ran[1] * diff.y + ll.y;
			point.z = ran[2] * diff.z + ll.z;
			res.push_back(point);
		}
	}
	else
	{
		int rank = 0;
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		for (size_t i = 0; i < PointNum; ++i)
		{
			ran[0] = dist(generator);
			ran[1] = dist(generator);
			ran[2] = dist(generator);
			point.x = ran[0] * diff.x + ll.x;
			point.y = ran[1] * diff.y + ll.y;
			point.z = ran[2] * diff.z + ll.z;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			  if (PointInPoly(*tproc, point, static_cast<size_t>(rank)))
					res.push_back(point);
		}
	}
	return res;
}

vector<Vector3D> RandRectangular(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, boost::mt19937_64 &generator)
{
	double ran[3];
	Vector3D diff = ur - ll;
	vector<Vector3D> res;
	Vector3D point;
	boost::random::uniform_real_distribution<> dist;
	for (size_t i = 0; i < PointNum; ++i)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		ran[2] = dist(generator);
		point.x = ran[0] * diff.x + ll.x;
		point.y = ran[1] * diff.y + ll.y;
		point.z = ran[2] * diff.z + ll.z;
		res.push_back(point);
	}
	return res;
}


vector<Vector3D> RandSphereR2(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax, const Vector3D& center,
	Voronoi3D const* tproc)
{
	typedef boost::mt19937_64 base_generator_type;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	vector<Vector3D> res;
	if (tproc == nullptr)
	{
		res.reserve(PointNum);
		while (res.size() < PointNum)
		{
			double r = dist(generator)*(Rmax - Rmin) + Rmin;
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
				res.push_back(point);
		}
	}
	else
	{
		int rank = 0;
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		for (size_t i = 0; i < PointNum; ++i)
		{
			double r = dist(generator)*(Rmax - Rmin) + Rmin;
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			  if (PointInPoly(*tproc, point, static_cast<size_t>(rank)))
					res.push_back(point);
		}
	}
	return res;
}

vector<Vector3D> RandSphereR(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax,
	const Vector3D& center, Voronoi3D const* tproc)
{
	typedef boost::mt19937_64 base_generator_type;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	vector<Vector3D> res;
	if (tproc == nullptr)
	{
		res.reserve(PointNum);
		while (res.size() < PointNum)
		{
			double r = std::pow(dist(generator)*(Rmax*Rmax*Rmax - Rmin * Rmin*Rmin) + Rmin * Rmin*Rmin, 0.333333333);
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
				res.push_back(point);
		}
	}
	else
	{
		int rank = 0;
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		for (size_t i = 0; i < PointNum; ++i)
		{
			double r = std::pow(dist(generator)*(Rmax*Rmax*Rmax - Rmin * Rmin*Rmin) + Rmin * Rmin*Rmin, 0.333333333);
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			  if (PointInPoly(*tproc, point, static_cast<size_t>(rank)))
					res.push_back(point);
		}
	}
	return res;
}

vector<Vector3D> RandSphereR1(std::size_t PointNum, Vector3D const& ll, Vector3D const& ur, double Rmin, double Rmax, const Vector3D& center,
	Voronoi3D const* tproc)
{
	typedef boost::mt19937_64 base_generator_type;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	vector<Vector3D> res;
	if (tproc == nullptr)
	{
		res.reserve(PointNum);
		while (res.size() < PointNum)
		{
			double r = std::sqrt(dist(generator)*(Rmax - Rmin)*(Rmax - Rmin)) + Rmin;
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
				res.push_back(point);
		}
	}
	else
	{
		int rank = 0;
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		for (size_t i = 0; i < PointNum; ++i)
		{
			double r = std::sqrt(dist(generator)*(Rmax - Rmin)*(Rmax - Rmin)) + Rmin;
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			  if (PointInPoly(*tproc, point, static_cast<size_t>(rank)))
					res.push_back(point);
		}
	}
	return res;
}

vector<Vector3D> RandSphereRa(std::size_t PointNum, Vector3D const & ll, Vector3D const & ur, double Rmin, double Rmax, double a, Vector3D const& center,
	Voronoi3D const* tproc)
{
	typedef boost::mt19937_64 base_generator_type;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	vector<Vector3D> res;
	double Rmx = std::pow(Rmax, a);
	double Rmn = std::pow(Rmin, a);
	double a_1 = 1.0 / a;
	if (tproc == nullptr)
	{
		res.reserve(PointNum);
		while (res.size() < PointNum)
		{
			double r = std::pow(dist(generator)*(Rmx - Rmn) + Rmn, a_1);
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
				res.push_back(point);
		}
	}
	else
	{
		int rank = 0;
#ifdef RICH_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
		for (size_t i = 0; i < PointNum; ++i)
		{
			double r = std::pow(dist(generator)*(Rmx - Rmn) + Rmn, a_1);
			double phi = 2 * M_PI*dist(generator);
			double t = acos(2 * dist(generator) - 1);
			Vector3D point(r*sin(t)*cos(phi), r*sin(t)*sin(phi), r*cos(t));
			point += center;
			if (point.x<ur.x&&point.x>ll.x&&point.y > ll.y&&point.y<ur.y&&point.z>ll.z&&point.z < ur.z)
			  if (PointInPoly(*tproc, point, static_cast<size_t>(rank)))
					res.push_back(point);
		}
	}
	return res;
}


#ifdef RICH_MPI
vector<Vector3D> RandPointsMPI(Voronoi3D const& tproc, size_t np)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	boost::mt19937_64 generator(rank);
	face_vec const& faces = tproc.GetCellFaces(rank);
	point_vec points = tproc.GetPointsInFace(faces[0]);
	Vector3D ll = tproc.GetFacePoints()[points[0]];
	Vector3D ur(ll);
	for (size_t i = 0; i < faces.size(); ++i)
	{
		points = tproc.GetPointsInFace(faces[i]);
		for (size_t j = 0; j < points.size(); ++j)
		{
			ll.x = std::min(ll.x, tproc.GetFacePoints()[points[j]].x);
			ll.y = std::min(ll.y, tproc.GetFacePoints()[points[j]].y);
			ll.z = std::min(ll.z, tproc.GetFacePoints()[points[j]].z);
			ur.x = std::max(ur.x, tproc.GetFacePoints()[points[j]].x);
			ur.y = std::max(ur.y, tproc.GetFacePoints()[points[j]].y);
			ur.z = std::max(ur.z, tproc.GetFacePoints()[points[j]].z);
		}
	}
	double v = tproc.GetVolume(rank);
	double TotVolume = (tproc.GetBoxCoordinates().second.x - tproc.GetBoxCoordinates().first.x)*
		(tproc.GetBoxCoordinates().second.y - tproc.GetBoxCoordinates().first.y)*
		(tproc.GetBoxCoordinates().second.z - tproc.GetBoxCoordinates().first.z);
	np = static_cast<size_t>(floor(v*static_cast<double>(np) / TotVolume + 0.5));
	vector<Vector3D> res;
	double ran[3];
	boost::random::uniform_real_distribution<> dist;
	while (res.size() < np)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		ran[2] = dist(generator);
		Vector3D test_point(ll.x + ran[0] * (ur.x - ll.x), ll.y + ran[1] * (ur.y - ll.y), ll.z + ran[2] * (ur.z - ll.z));
		if (PointInPoly(tproc, test_point, rank))
			res.push_back(test_point);
	}
	return res;
}
#endif


