/*! \file Tessellation3D.hpp
\brief Class for the Cube cartesian tessellation in 3D
\author Elad Steinberg
*/

#ifndef CUBE3D_HPP
#define CUBE3D_HPP 1

#include "Tessellation3D.hpp"
#include "OuterBoundary3D.hpp"
#include "../../misc/utils.hpp"

class CubeGrid : public Tessellation3D
{
private:
	size_t nx_,ny_,nz_;
	double dx_,dy_,dz_;
	OuterBoundary3D const* obc_;
	Vector3D backlowerleft_,frontupperright_;
	vector<Vector3D> cor_;
	vector<Face> faces_;
	vector<vector<size_t> > cellfaces_;
public:
	CubeGrid(void);

	CubeGrid(CubeGrid const& other);

	void Initialise(vector<Vector3D> const& points, OuterBoundary3D const* bc);

	void Update(vector<Vector3D> const& points);

	size_t GetPointNo(void) const;

	Vector3D GetMeshPoint(size_t index) const;

	Vector3D const& GetCellCM(size_t index) const;

	size_t GetTotalFacesNumber(void) const;

	Face const& GetEdge(size_t index) const;

	double GetWidth(size_t index) const;

	double GetVolume(size_t index) const;

	vector<size_t>const& GetCellFaces(int index) const;

	vector<Vector3D>& GetMeshPoints(void);

	vector<size_t> GetNeighbors(size_t index)const;

	CubeGrid* clone(void) const;

	~CubeGrid(void);

	bool NearBoundary(size_t index) const;

	vector<vector<size_t> >& GetDuplicatedPoints(void);

	vector<vector<size_t> >const& GetDuplicatedPoints(void)const;

	size_t GetTotalPointNumber(void)const;

	vector<Vector3D>& GetAllCM(void);

	void GetNeighborNeighbors(vector<size_t> &result, int point)const;
};


#endif //CUBE3D_HPP