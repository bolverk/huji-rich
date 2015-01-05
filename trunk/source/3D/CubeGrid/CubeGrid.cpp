#include "CubeGrid.hpp"

CubeGrid::CubeGrid(void):nx_(0),ny_(0),nz_(0),dx_(0),dy_(0),dz_(0),obc_(0),
	backlowerleft_(Vector3D()),frontupperright_(Vector3D()),cor_(vector<Vector3D> ()),
	faces_(vector<Face> ()),cellfaces_(vector<vector<size_t> > ()){}

CubeGrid::CubeGrid(CubeGrid const& other):nx_(other.nx_),ny_(other.ny_),nz_(other.nz_),
	dx_(other.dx_),dy_(other.dy_),dz_(other.dz_),obc_(other.obc_),
	backlowerleft_(other.backlowerleft_),frontupperright_(other.frontupperright_),
	cor_(other.cor_),faces_(other.faces_),cellfaces_(other.cellfaces_){}

void CubeGrid::Initialise(vector<Vector3D> const& points, OuterBoundary3D const* bc)
{
	obc_=bc;
	backlowerleft_=obc_->GetGridBoundary(BackLowerLeft);
	frontupperright_=obc_->GetGridBoundary(FrontUpperRight);

/*	const size_t nz=index%(nz_);
	const size_t nx=index/(nz_*ny_);
	const size_t ny=(index-nx*nz_*ny_)/nz_;
	return Vector3D(backlowerleft_.x+(nx_+0.5)*dx_,backlowerleft_.y+(ny_+0.5)*dy_,
		backlowerleft_.z+(nz_+0.5)*dz_);
		*/
	cor_=points;
	// Create the faces
	size_t nfaces=GetTotalFacesNumber();
	faces_.resize(nfaces);
	cellfaces_.resize(cor_.size());
	for(size_t i=0;i<cellfaces_.size();++i)
		cellfaces_[i].resize(6);

	// Do the z direction
	//	size_t ntemp=(nx_+1)*ny_+(ny_+1)*nx_;
	vector<Vector3D> side(4);
	side[0]=backlowerleft_;
	side[1]=backlowerleft_+Vector3D(dx_,0,0);
	side[2]=backlowerleft_+Vector3D(dx_,dy_,0);
	side[3]=backlowerleft_+Vector3D(0,dy_,0);

	for(size_t i=0;i<ny_+1;++i)
	for(size_t j=0;j<nx_+1;++j)
	{
		faces_[i*(nx_+1)+j].vertices=side+Vector3D(dx_*j,0,0)+Vector3D(0,dy_*i,0);
	}
}

void CubeGrid::Update(vector<Vector3D> const& /*points*/)
{}

size_t CubeGrid::GetPointNo(void) const
{
	return nz_*nx_*ny_;
}

Vector3D CubeGrid::GetMeshPoint(size_t index) const
{
	return cor_[index];
}

Vector3D const& CubeGrid::GetCellCM(size_t index) const
{
  return cor_[index];
}

size_t CubeGrid::GetTotalFacesNumber(void) const
{
	return (nz_+1)*nx_*ny_+(nx_+1)*nz_*ny_+(ny_+1)*nx_*nz_;
}

double CubeGrid::GetWidth(size_t /*index*/) const
{
	return min(min(dx_,dy_),dz_);
}

double CubeGrid::GetVolume(size_t /*index*/) const
{
	return dx_*dy_*dz_;
}

vector<size_t>const& CubeGrid::GetCellFaces(int /*index*/) const
{
  throw;
}
