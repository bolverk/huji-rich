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
#include <boost/array.hpp>
#include "tetgen.h"
#include "Vector3D.hpp"
#include <stack>
#include <set>
#include "Tessellation3D.hpp"
#include "Intersections.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

typedef boost::array<int, 4> b_array_4;
typedef boost::array<int, 3> b_array_3;

class Tetrahedron
{
public:
	b_array_4 points, neighbors;

	Tetrahedron(void);
};


class Voronoi3D : public Tessellation3D
{
private:
	Vector3D ll_, ur_;
	size_t Norg_, bigtet_;

	std::set<int> set_temp_;
	std::stack<int> stack_temp_;

#ifdef RICH_MPI
	vector<Vector3D> UpdateMPIPoints(Tessellation3D const& vproc, int rank,
		vector<Vector3D> const& points, vector<size_t> &selfindex, vector<int> &sentproc,
		vector<vector<size_t> > &sentpoints);
#endif
	vector<size_t> FindIntersectionsRecursive(Tessellation3D const& tproc, size_t rank, Sphere const& sphere);
	size_t GetFirstPointToCheck(void)const;
	void GetPointToCheck(size_t point, vector<bool> const& checked, vector<size_t> &res);
	void CalcRigidCM(size_t face_index);
	void RunTetGen(vector<Vector3D> const& points,tetgenio &tetin,tetgenio &tetout, bool voronoi = false);
	Vector3D GetTetraCM(boost::array<Vector3D, 4> const& points)const;
	double GetTetraVolume(boost::array<Vector3D, 4> const& points)const;
	void CalcCellCMVolume(size_t index);
	double GetRadius(size_t index);
	double GetMaxRadius(size_t index);
	vector<std::pair<size_t, size_t> > FindIntersections(
#ifdef RICH_MPI
		Tessellation3D const& tproc
#endif
		);
	void CopyData(tetgenio &tetin);
	void CopyDataVoronoi(tetgenio &tetin);
	void FillPointTetra(size_t point, size_t initetra);
	double CalcTetraRadiusCenter(size_t index);
	vector<Vector3D> CreateBoundaryPoints(vector<std::pair<size_t, size_t> > const& to_duplicate);
	vector<Vector3D> CreateBoundaryPointsMPI(vector<std::pair<size_t, size_t> > const& to_duplicate,
		Tessellation3D const& tproc);
	vector<Vector3D> mesh_points_;
	vector<Tetrahedron> tetras_;
	vector<vector<size_t> > PointTetras_; // The tetras containing each point
	vector<double> R_; // The radius of the sphere of each tetra
	vector<Vector3D> tetra_centers_;
	// Voronoi Data
	vector<vector<size_t> > FacesInCell_;
	vector<vector<size_t> > PointsInFace_;
	vector<std::pair<size_t, size_t> > FaceNeighbors_;
	vector<Vector3D> CM_;
	vector<double> volume_;
	vector<double> area_;
	vector<Vector3D> FacePoints_;
	vector<vector<size_t> > duplicated_points_;
	vector<int> sentprocs_,duplicatedprocs_;
	vector<vector<size_t> > sentpoints_;
	vector<size_t> self_index_,self_duplicate_;
	Voronoi3D();
public:
	Vector3D FaceCM(size_t index)const;

	Voronoi3D(Vector3D const& ll, Vector3D const& ur);

	void output(std::string const& filename)const;

	void Build(vector<Vector3D> const& points
#ifdef RICH_MPI
		, Tessellation3D const& tproc
#endif
		);

	size_t GetPointNo(void) const;

	Vector3D GetMeshPoint(size_t index) const;

	double GetArea(size_t index) const;

	Vector3D const& GetCellCM(size_t index) const;

	size_t GetTotalFacesNumber(void) const;

	double GetWidth(size_t index) const;

	double GetVolume(size_t index) const;

	vector<size_t>const& GetCellFaces(size_t index) const;
	
	vector<Vector3D>& GetMeshPoints(void);

	vector<size_t> GetNeighbors(size_t index)const;

	Tessellation3D* clone(void) const;

	bool NearBoundary(size_t index) const;

	bool BoundaryFace(size_t index) const;

	vector<vector<size_t> >& GetDuplicatedPoints(void);

	vector<vector<size_t> >const& GetDuplicatedPoints(void)const;

	size_t GetTotalPointNumber(void)const;

	vector<Vector3D>& GetAllCM(void);

	void GetNeighborNeighbors(vector<size_t> &result, size_t point)const;

	Vector3D Normal(size_t faceindex)const;

	bool IsGhostPoint(size_t index)const;

	Vector3D CalcFaceVelocity(size_t index, Vector3D const& v0, Vector3D const& v1)const;

	vector<Vector3D>const& GetFacePoints(void) const;

	vector<size_t>const& GetPointsInFace(size_t index) const;

	std::pair<size_t,size_t> GetFaceNeighbors(size_t face_index)const;

	vector<int> GetDuplicatedProcs(void)const;

	vector<int> GetSentProcs(void)const;

	vector<vector<size_t> > const& GetSentPoints(void)const;

	vector<size_t> const& GetSelfIndex(void) const;

	vector<size_t> const& GetSelfDuplicate(void)const;
};

#endif // VORONOI3D_HPP