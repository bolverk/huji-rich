/*! \file Tessellation3D.hpp
\brief Class for the Cube cartesian tessellation in 3D
\author Elad Steinberg
*/

#ifndef CUBE3D_HPP
#define CUBE3D_HPP 1

#include "../GeometryCommon/Tessellation3D.hpp"
#include "../GeometryCommon/OuterBoundary3D.hpp"
#include "../../misc/utils.hpp"

//! \brief Cartesian 3d grid
class CubeGrid : public Tessellation3D
{
private:
	const size_t nx_,ny_,nz_,maxsize_;
	const double dx_,dy_,dz_;
	OuterBoundary3D const* obc_;
	Vector3D backlowerleft_,frontupperright_;
	vector<Vector3D> cor_;
	vector<Face> faces_;
	vector<vector<size_t> > cellfaces_;
	vector<vector<size_t> > temp_;

	CubeGrid(CubeGrid const& other);
	CubeGrid& operator=(CubeGrid const& other);
public:
	/*!
	\brief Class constructor
	\param nx Number of points in the x direction
	\param ny Number of points in the y direction
	\param nz Number of points in the z direction
	\param backlowerleft The boundary of the domain with the lowest x,y,z
	\param frontupperright The boundary of the domain with the highest x,y,z
	*/
	CubeGrid(size_t nx, size_t ny, size_t nz, Vector3D const& backlowerleft,
		Vector3D const& frontupperright);

	void Initialise(vector<Vector3D> const& points, OuterBoundary3D const* bc);

	void Update(vector<Vector3D> const& points);

	size_t GetPointNo(void) const;

	Vector3D GetMeshPoint(size_t index) const;

	Vector3D const& GetCellCM(size_t index) const;

	size_t GetTotalFacesNumber(void) const;

	Face const& GetFace(size_t index) const;

	double GetWidth(size_t index) const;

	double GetVolume(size_t index) const;

	vector<size_t>const& GetCellFaces(size_t index) const;

	vector<Vector3D>& GetMeshPoints(void);

	vector<size_t> GetNeighbors(size_t index)const;

	Tessellation3D* clone(void) const;

	~CubeGrid(void);

	bool NearBoundary(size_t index) const;

	vector<vector<size_t> >& GetDuplicatedPoints(void);

	vector<vector<size_t> >const& GetDuplicatedPoints(void)const;

	size_t GetTotalPointNumber(void)const;

	vector<Vector3D>& GetAllCM(void);

	void GetNeighborNeighbors(vector<size_t> &result, size_t point)const;

	Vector3D Normal(size_t faceindex)const;

	bool IsGhostPoint(size_t index)const;

	Vector3D CalcFaceVelocity(size_t p0,size_t p1,Vector3D const& v0,
		Vector3D const& v1)const;

	bool BoundaryFace(size_t index) const;
};


#endif //CUBE3D_HPP
