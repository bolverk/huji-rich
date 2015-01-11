#include "CubeGrid.hpp"
#include <limits>

namespace
{
	/*void ConvertIndexToSubindex(size_t index,size_t &i,size_t &j,size_t &k,
		size_t ny_,size_t nz_)
	{
		k=index%(nz_);
		i=index/(nz_*ny_);
		j=(index-i*nz_*ny_)/nz_;
	}
	*/
	size_t ConvertSubindexToIndex(size_t i,size_t j,size_t k,
		size_t ny_,size_t nz_)
	{
		return i*ny_*nz_+j*nz_+k;
	}
}

CubeGrid::CubeGrid(size_t nx, size_t ny, size_t nz, Vector3D const& backlowerleft,
	Vector3D const& frontupperright):nx_(nx),ny_(ny),nz_(nz),
	maxsize_(std::numeric_limits<std::size_t>::max()),
	dx_((frontupperright.x-backlowerleft.x)/nx),
	dy_((frontupperright.y-backlowerleft.y)/ny),
	dz_((frontupperright.z-backlowerleft.z)/nz),
	obc_(0),backlowerleft_(backlowerleft),frontupperright_(frontupperright),
	cor_(vector<Vector3D> ()),faces_(vector<Face> ()),cellfaces_(vector<vector<size_t> > ()),
	temp_(vector<vector<size_t> > ())
{
	// Create the faces
	size_t nfaces=(nx_+1)*ny_*nz_+(ny_+1)*nx_*nz_+(nz_+1)*ny_*nx_;
	faces_.resize(nfaces);
	cellfaces_.resize(nx_*ny_*nz_);

	// Do the z direction
	vector<Vector3D> side(4);
	side[0]=backlowerleft_;
	side[1]=backlowerleft_+Vector3D(dx_,0,0);
	side[2]=backlowerleft_+Vector3D(dx_,dy_,0);
	side[3]=backlowerleft_+Vector3D(0,dy_,0);

	for(size_t k=0;k<nz_+1;++k)
	{
		for(size_t i=0;i<ny_;++i)
		{
			for(size_t j=0;j<nx_;++j)
			{
				faces_[i*nx_+j+k*ny_*nx_].vertices=side+Vector3D(dx_*j,0,dz_*k)
					+Vector3D(0,dy_*i,0);
				if(k==0)
					faces_[i*nx_+j+k*ny_*nx_].neighbors.first=maxsize_;
				else
					faces_[i*nx_+j+k*ny_*nx_].neighbors.first=ConvertSubindexToIndex(j,i,k-1,ny_,nz_);
				if(k==nz_)
					faces_[i*nx_+j+k*ny_*nx_].neighbors.second=maxsize_;
				else
					faces_[i*nx_+j+k*ny_*nx_].neighbors.second=ConvertSubindexToIndex(j,i,k,ny_,nz_);
			}
		}
	}

	// Do the y direction
	side[0]=backlowerleft_;
	side[1]=backlowerleft_+Vector3D(dx_,0,0);
	side[2]=backlowerleft_+Vector3D(dx_,0,dz_);
	side[3]=backlowerleft_+Vector3D(0,0,dz_);

	for(size_t k=0;k<nz_;++k)
	{
		for(size_t i=0;i<ny_+1;++i)
		{
			for(size_t j=0;j<nx_;++j)
			{
				faces_[i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_].vertices=side+Vector3D(dx_*j,dy_*i,0)
					+Vector3D(0,0,dz_*k);
				if(i==0)
					faces_[i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_].neighbors.first=maxsize_;
				else
					faces_[i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_].neighbors.first=ConvertSubindexToIndex(j,i-1,k,ny_,nz_);
				if(i==ny_)
					faces_[i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_].neighbors.second=maxsize_;
				else
					faces_[i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_].neighbors.second=ConvertSubindexToIndex(j,i,k,ny_,nz_);
			}
		}
	}
	// Do the x direction
	side[0]=backlowerleft_;
	side[1]=backlowerleft_+Vector3D(0,dy_,0);
	side[2]=backlowerleft_+Vector3D(0,dy_,dz_);
	side[3]=backlowerleft_+Vector3D(0,0,dz_);

	for(size_t k=0;k<nz_;++k)
	{
		for(size_t i=0;i<ny_;++i)
		{
			for(size_t j=0;j<nx_+1;++j)
			{
				faces_[i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_].vertices=side
					+Vector3D(0,dy_*i,0)+Vector3D(dx_*j,0,dz_*k);
				if(j==0)
					faces_[i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_].neighbors.first=maxsize_;
				else
					faces_[i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_].neighbors.first=ConvertSubindexToIndex(j-1,i,k,ny_,nz_);
				if(j==nx_)
					faces_[i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_].neighbors.second=maxsize_;
				else
					faces_[i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_].neighbors.second=ConvertSubindexToIndex(j,i,k,ny_,nz_);
			}
		}
	}
	// Add the faces for the points
	for(size_t i=0;i<cellfaces_.size();++i)
		cellfaces_[i].reserve(6);
	for(size_t k=0;k<nz_;++k)
	{
		for(size_t i=0;i<ny_;++i)
		{
			for(size_t j=0;j<nx_;++j)
			{
				size_t index=ConvertSubindexToIndex(j,i,k,ny_,nz_);
				cellfaces_[index].push_back(i*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_+nx_*(nz_+1)*ny_);
				cellfaces_[index].push_back((i+1)*nx_*nz_+k+j*nz_+nx_*(nz_+1)*ny_+nx_*(nz_+1)*ny_);

				cellfaces_[index].push_back(i*nz_+j*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_);
				cellfaces_[index].push_back(i*nz_+(j+1)*ny_*nz_+k+nx_*(nz_+1)*ny_+nx_*(ny_+1)*nz_);

				cellfaces_[index].push_back(i*nx_+j+k*ny_*nx_);
				cellfaces_[index].push_back(i*nx_+j+(k+1)*ny_*nx_);
			}
		}
	}
	// Create the points
	cor_.resize(nx_*ny_*nz_);
	for(size_t i=0;i<nx_;++i)
	{
		for(size_t j=0;j<ny_;++j)
		{
			for(size_t k=0;k<nz_;++k)
			{
				cor_[ConvertSubindexToIndex(i,j,k,ny_,nz_)]=backlowerleft+
					Vector3D(i*dx_+0.5*dx_,j*dy_+0.5*dy_,k*dz_+0.5*dz_);
			}
		}
	}
}

CubeGrid::CubeGrid(CubeGrid const& other):nx_(other.nx_),ny_(other.ny_),nz_(other.nz_),
	maxsize_(other.maxsize_),dx_(other.dx_),dy_(other.dy_),dz_(other.dz_),obc_(other.obc_),
	backlowerleft_(other.backlowerleft_),frontupperright_(other.frontupperright_),
	cor_(other.cor_),faces_(other.faces_),cellfaces_(other.cellfaces_),
	temp_(other.temp_){}

void CubeGrid::Initialise(vector<Vector3D> const& /*points*/, OuterBoundary3D const* bc)
{
	obc_=bc;
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
	return faces_.size();
}

Face const& CubeGrid::GetFace(size_t index) const
{
	return faces_[index];
}

double CubeGrid::GetWidth(size_t /*index*/) const
{
	return min(min(dx_,dy_),dz_);
}

double CubeGrid::GetVolume(size_t /*index*/) const
{
	return dx_*dy_*dz_;
}

vector<size_t>const& CubeGrid::GetCellFaces(size_t index) const
{
	return cellfaces_[index];
}

vector<Vector3D>& CubeGrid::GetMeshPoints(void)
{
	return cor_;
}

vector<size_t> CubeGrid::GetNeighbors(size_t index)const
{
	vector<size_t> res(6);
	for(size_t i=0;i<6;++i)
	{
		res[i]=(faces_[cellfaces_[index][i]].neighbors.first==index) ?
			faces_[cellfaces_[index][i]].neighbors.second :
			faces_[cellfaces_[index][i]].neighbors.first;
	}
	return res;
}

Tessellation3D* CubeGrid::clone(void) const
{
	return new CubeGrid(*this);
}

CubeGrid::~CubeGrid(void)
{}

bool CubeGrid::NearBoundary(size_t index) const
{
	for(size_t i=0;i<6;++i)
	{
		if(faces_[cellfaces_[index][i]].neighbors.first==maxsize_)
			return true;
		if(faces_[cellfaces_[index][i]].neighbors.second==maxsize_)
			return true;
	}
	return false;
}

vector<vector<size_t> >& CubeGrid::GetDuplicatedPoints(void)
{
	return temp_;
}

vector<vector<size_t> >const& CubeGrid::GetDuplicatedPoints(void)const
{
	return temp_;
}

size_t CubeGrid::GetTotalPointNumber(void)const
{
	return cor_.size();
}

vector<Vector3D>& CubeGrid::GetAllCM(void)
{
	return cor_;
}

void CubeGrid::GetNeighborNeighbors(vector<size_t> &result, size_t point)const
{
	result.clear();
	vector<size_t> neigh=GetNeighbors(point);
	for(size_t i=0;i<neigh.size();++i)
	{
		if(neigh[i]==maxsize_)
			continue;
		vector<size_t> temp=GetNeighbors(neigh[i]);
		result.insert(result.end(),temp.begin(),temp.end());
	}
	sort(result.begin(),result.end());
	result=unique(result);
}

Vector3D CubeGrid::Normal(size_t faceindex)const
{
	if(faceindex<(nx_*ny_*(nz_+1)))
		return Vector3D(0,0,dz_);
	if(faceindex<(nx_*(ny_+1)*(nz_+1)))
		return Vector3D(0,dy_,0);
	else
		return Vector3D(dx_,0,0);
}

bool CubeGrid::IsGhostPoint(size_t index)const
{
	if(index<nz_*ny_*nx_)
		return false;
	else
		return true;
}

Vector3D CubeGrid::CalcFaceVelocity(size_t /*p0*/,size_t /*p1*/,Vector3D const& /*v0*/,
		Vector3D const& /*v1*/)const
{
	return Vector3D();
}

bool CubeGrid::BoundaryFace(size_t index) const
{
	if(faces_[index].neighbors.first==maxsize_||faces_[index].neighbors.second==maxsize_)
		return true;
	else
		return false;
}
