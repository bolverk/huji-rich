#include "RoundGrid.hpp"
#include "HilbertOrder.hpp"

namespace {
  Vector2D genNewPoint
  (const Tessellation& tess,
   int i)
  {
    const double eta_=0.02,chi_=1;
    const double R = sqrt(tess.GetVolume(i)/M_PI);
    const Vector2D s = tess.GetCellCM(i);
    const Vector2D r = tess.GetMeshPoint(i);
    const double d = abs(s-r);
    const Vector2D dw =
      d/eta_/R<0.95 ? Vector2D(0,0) : chi_*0.5*(s-r);
    return tess.GetMeshPoint(i) + dw;
  }
}

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
  vector<Vector2D> res
    (static_cast<size_t>(tess->GetPointNo()));
  for(int j=0,N=tess->GetPointNo();j<NumberIt;++j)
    {
#ifdef RICH_MPI
      N=tess->GetPointNo();
      res.resize(static_cast<size_t>(N));
#endif
      for(int i=0;i<N;++i)
	  res[static_cast<size_t>(i)] = genNewPoint(*tess, i);
#ifdef RICH_MPI
      tproc == 0 ? tess->Update(res) : tess->Update(res,*tproc);
#else
      tess->Update(res);
#endif
    }
#ifdef RICH_MPI
  res=tess->GetMeshPoints();
  res.resize(static_cast<size_t>(tess->GetPointNo()));
#endif
  return res;
}
