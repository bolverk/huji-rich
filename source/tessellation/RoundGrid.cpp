#include "RoundGrid.hpp"
#include "HilbertOrder.hpp"

vector<Vector2D> RoundGrid(vector<Vector2D> const& points,
			   const OuterBoundary& bc,int NumberIt,
			   #ifdef RICH_MPI
			   Tessellation const* tproc,
			   #endif
	Tessellation *tess)
{
  VoronoiMesh default_tess;
  tess = tess == nullptr ? &default_tess : tess;
#ifdef RICH_MPI
  tproc == 0 ? tess->Initialise(points,bc) : tess->Initialise(points,*tproc,bc);
#else
	tess->Initialise(points,bc);
#endif
	const double eta_=0.02,chi_=1;
	vector<Vector2D> res
	  (static_cast<size_t>(tess->GetPointNo()));
	for(int j=0,N=tess->GetPointNo();j<NumberIt;++j)
	{
#ifdef RICH_MPI
	  N=tess->GetPointNo();
		res=tess->GetMeshPoints();
		res.resize(static_cast<size_t>(N));
#endif
		for(int i=0;i<N;++i)
		{
			double R = sqrt(tess->GetVolume(i)/M_PI);
			Vector2D s = tess->GetCellCM(i);
			Vector2D r = tess->GetMeshPoint(i);
			double d = abs(s-r);
			Vector2D dw;
			if(d/eta_/R<0.95)
				dw = 0*s;
			else
				dw = chi_*0.5*(s-r);
			res[static_cast<size_t>(i)]=tess->GetMeshPoint(i)+dw;
		}
#ifdef RICH_MPI
		tproc == 0 ? tess->Update(res) : tess->Update(res,*tproc);
#else
		tess->Update(res);
#endif
	}
#ifdef RICH_MPI
	const int N=tess->GetPointNo();
	res=tess->GetMeshPoints();
	res.resize(static_cast<size_t>(N));
#endif
	return res;
}
