/* \file Voronoi3D.hpp
\brief A 3D Voronoi
\Author Elad Steinberg
*/
#ifndef VORONOI3D_HPP
#define VORONOI3D_HPP 1

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include "Delaunay3D.hpp"
#include "Intersections.hpp"
#include <stack>
#include <set>
#include <boost/array.hpp>
#include "Tessellation3D.hpp"
#include <boost/container/flat_set.hpp>

#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

typedef boost::array<std::size_t, 4> b_array_4;
typedef boost::array<std::size_t, 3> b_array_3;

class Voronoi3D : public Tessellation3D
{
private:
	Vector3D ll_, ur_;
	std::size_t Norg_, bigtet_;

	std::set<int> set_temp_;
	std::stack<int> stack_temp_;

#ifdef RICH_MPI
	vector<Vector3D> UpdateMPIPoints(Tessellation3D const& vproc, int rank,
		vector<Vector3D> const& points, vector<std::size_t> &selfindex, vector<int> &sentproc,
		vector<vector<std::size_t> > &sentpoints);
#endif
	void FindIntersectionsSingle(vector<Face> const& box, std::size_t point, Sphere &sphere,
		vector<size_t> &intersecting_faces);
	vector<std::size_t> FindIntersectionsRecursive(Tessellation3D const& tproc, std::size_t rank, 
		std::size_t point, Sphere &sphere, bool recursive, boost::container::flat_set<size_t> &visited,
		std::stack<std::size_t> &to_check);
	std::size_t GetFirstPointToCheck(void)const;
	void GetPointToCheck(std::size_t point, vector<bool> const& checked, vector<std::size_t> &res);
	void CalcRigidCM(std::size_t face_index);
	Vector3D GetTetraCM(boost::array<Vector3D, 4> const& points)const;
	double GetTetraVolume(boost::array<Vector3D, 4> const& points)const;
	void CalcCellCMVolume(std::size_t index);
	double GetRadius(std::size_t index);
	double GetMaxRadius(std::size_t index);
	void CalcAllCM(void);
	vector<std::pair<std::size_t, std::size_t> > SerialFindIntersections(bool first_run);
	vector<std::pair<std::size_t, std::size_t> > SerialFirstIntersections(void);
#ifdef RICH_MPI
	vector<std::pair<std::size_t, std::size_t> > FindIntersections(Tessellation3D const& tproc, bool recursive);
	vector<Vector3D> CreateBoundaryPointsMPI(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
		Tessellation3D const& tproc, vector<vector<size_t> > &self_duplicate);
	//void CleanDuplicate(vector<std::pair<size_t, size_t> > const& ghost_index, Tessellation3D const& tproc);
#endif
	double CalcTetraRadiusCenter(std::size_t index);
	vector<Vector3D> CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
		vector<vector<size_t> > &past_duplicate);
	void BuildVoronoi(void);

	Delaunay3D del_;
	vector<vector<std::size_t> > PointTetras_; // The tetras containing each point
	vector<double> R_; // The radius of the sphere of each tetra
	vector<Vector3D> tetra_centers_;
	// Voronoi Data
	vector<vector<std::size_t> > FacesInCell_;
	vector<vector<std::size_t> > PointsInFace_; // Right hand with regard to first neighbor
	vector<std::pair<std::size_t, std::size_t> > FaceNeighbors_;
	vector<Vector3D> CM_,Face_CM_;
	vector<double> volume_;
	vector<double> area_;
	vector<vector<std::size_t> > duplicated_points_;
	vector<int> sentprocs_, duplicatedprocs_;
	vector<vector<std::size_t> > sentpoints_, Nghost_;
	vector<std::size_t> self_index_;
	Voronoi3D();
public:
	vector<Vector3D>& GetAllFaceCM(void);

	Vector3D FaceCM(std::size_t index)const;

	Voronoi3D(Vector3D const& ll, Vector3D const& ur);

	void output(std::string const& filename)const;

	void Build(vector<Vector3D> const& points);

#ifdef RICH_MPI
	void Build(vector<Vector3D> const& points, Tessellation3D const& tproc);

	friend void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points, size_t Niter, double speed, int mode);
#endif

	std::size_t GetPointNo(void) const;

	Vector3D GetMeshPoint(std::size_t index) const;

	double GetArea(std::size_t index) const;

	Vector3D const& GetCellCM(std::size_t index) const;

	std::size_t GetTotalFacesNumber(void) const;

	double GetWidth(std::size_t index) const;

	double GetVolume(std::size_t index) const;

	vector<std::size_t>const& GetCellFaces(std::size_t index) const;

	vector<Vector3D>& GetMeshPoints(void);

	vector<std::size_t> GetNeighbors(std::size_t index)const;

	Tessellation3D* clone(void) const;

	bool NearBoundary(std::size_t index) const;

	bool BoundaryFace(std::size_t index) const;

	vector<vector<std::size_t> >& GetDuplicatedPoints(void);

	vector<vector<std::size_t> >const& GetDuplicatedPoints(void)const;

	std::size_t GetTotalPointNumber(void)const;

	vector<Vector3D>& GetAllCM(void);

	void GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point)const;

	Vector3D Normal(std::size_t faceindex)const;

	bool IsGhostPoint(std::size_t index)const;

	Vector3D CalcFaceVelocity(std::size_t index, Vector3D const& v0, Vector3D const& v1)const;

	vector<Vector3D>& GetFacePoints(void);

	vector<double>& GetAllArea(void);

	vector<Vector3D>const& GetFacePoints(void) const;

	vector<vector<size_t> >& GetAllCellFaces(void);

	vector<std::size_t>const& GetPointsInFace(std::size_t index) const;

	std::pair<std::size_t, std::size_t> GetFaceNeighbors(std::size_t face_index)const;

	vector<int> GetDuplicatedProcs(void)const;

	vector<int> GetSentProcs(void)const;

	vector<vector<std::size_t> > const& GetSentPoints(void)const;

	vector<std::size_t> const& GetSelfIndex(void) const;

	vector<vector<std::size_t> > const& GetGhostIndeces(void) const;

	vector<vector<std::size_t> >& GetGhostIndeces(void);

	void GetNeighbors(size_t index, vector<size_t> &res)const;

	std::pair<Vector3D, Vector3D> GetBoxCoordinates(void)const;

	void BuildNoBox(vector<Vector3D> const& points, vector<Vector3D> const& ghosts,vector<size_t> toduplicate);

	vector<double>& GetAllVolumes(void);

	std::vector<std::pair<size_t, size_t> >& GetAllFaceNeighbors(void);

	vector<vector<size_t> > & GetAllPointsInFace(void);

	size_t& GetPointNo(void);

	bool IsPointOutsideBox(size_t index)const;
};

bool PointInPoly(Tessellation3D const& tess, Vector3D const& point, std::size_t index);

#endif // VORONOI3D_HPP