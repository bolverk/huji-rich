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

TessellationHandler::~TessellationHandler(void) = default;

void SerialHandler::initialise
(Tessellation& tess,
 const vector<Vector2D>& points,
 const OuterBoundary& bc) const
{
  tess.Initialise(points, bc);
}

void SerialHandler::update
(Tessellation& tess, const vector<Vector2D>& points) const
{
  tess.Update(points);
}

#ifdef RICH_MPI

ParallelHandler::ParallelHandler(const Tessellation& meta):
  meta_(meta) {}

void ParallelHandler::initialise
(Tessellation& tess,
 const vector<Vector2D>& points,
 const OuterBoundary& bc) const
{
  tess.Initialise(points, meta_, bc);
}

void ParallelHandler::update
(Tessellation& tess,
 const vector<Vector2D>& points) const
{
  tess.Update(points, meta_);
}

#endif // RICH_MPI

vector<Vector2D> RoundGrid(const vector<Vector2D>& points,
			   const OuterBoundary& bc,
#ifdef RICH_MPI
			   const TessellationHandler& th,
#endif // RICH_MPI
			   int NumberIt)
{
  VoronoiMesh tess;
  return RoundGrid
    (points,
     bc,
     tess,
#ifdef RICH_MPI
     th,
#endif // RICH_MPI
     NumberIt);
}

vector<Vector2D> RoundGrid
(const vector<Vector2D>& points,
 const OuterBoundary& bc,
 Tessellation& tess,
#ifdef RICH_MPI
 const TessellationHandler& th,
#endif
 int NumberIt)
{
#ifndef RICH_MPI
  const SerialHandler th;
#endif //RICH_MPI
  th.initialise(tess, points, bc);
  for(int j=0;j<NumberIt;++j)
    th.update(tess, genNewPoints(tess));
  return genNewPoints(tess);
}
