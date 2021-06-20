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

  vector<Vector2D> genNewPoints(Tessellation& tess)
  {
    const int N = tess.GetPointNo();
    vector<Vector2D> res(static_cast<size_t>(N));
    for(int i=0;i<N;++i)
      res[static_cast<size_t>(i)] = genNewPoint(tess,i);
    return res;
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
    (static_cast<size_t>(tess->GetPointNo()));
  for(int j=0;j<NumberIt;++j)
    {
      const vector<Vector2D> res = genNewPoints(*tess);
#ifdef RICH_MPI
      tproc == 0 ? tess->Update(res) : tess->Update(res,*tproc);
#else
      tess->Update(res);
#endif
    }
  return genNewPoints(*tess);
}
