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

  class TessellationHandler
  {
  public:

    virtual void initialise
    (Tessellation& tess,
     const vector<Vector2D>& points,
     const OuterBoundary& bc) const = 0;

    virtual void update
    (Tessellation& tess,
     const vector<Vector2D>& points) const = 0;

    virtual ~TessellationHandler(void) = default;
  };

  class SerialHandler: public TessellationHandler
  {
  public:

    SerialHandler(){;}

    void initialise
    (Tessellation& tess,
     const vector<Vector2D>& points,
     const OuterBoundary& bc) const
    {
      tess.Initialise(points, bc);
    }

    void update
    (Tessellation& tess, 
     const vector<Vector2D>& points) const
    {
      tess.Update(points);  
    }
  };

#ifdef RICH_MPI
  class ParallelHandler: public TessellationHandler
  {
  public:

    ParallelHandler(const Tessellation& meta):
      meta_(meta) {}

    void initialise
    (Tessellation& tess,
     const vector<Vector2D>& points,
     const OuterBoundary& bc) const
    {
      tess.Initialise(points, meta_, bc);
    }

    void update
    (Tessellation& tess,
     const vector<Vector2D>& points) const
    {
      tess.Update(points, meta_);
    }

  private:
    const Tessellation& meta_;
  };
#endif // RICH_MPI

  // Generic function
  vector<Vector2D> RoundGridG
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
}

#ifdef RICH_MPI

vector<Vector2D> RoundGridV
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
 int NumberIt)
{
  VoronoiMesh tess;
  return RoundGridG
    (points,
     bc,
     tess,
     SerialHandler(),
     NumberIt);
}

#endif // RICH_MPI

vector<Vector2D> RoundGridV
(const vector<Vector2D>& points,
 const OuterBoundary& bc,
#ifdef RICH_MPI
 const Tessellation& meta,
#endif // RICH_MPI
 int NumberIt)
{
  VoronoiMesh tess;
  return RoundGridG
    (points,
     bc,
     tess,
#ifdef RICH_MPI
     ParallelHandler(meta),
#endif // RICH_MPI
     NumberIt);
}

vector<Vector2D> RoundGrid
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
 Tessellation& tess,
#ifdef RICH_MPI
 const Tessellation& meta,
#endif
 int NumberIt)
{
  return RoundGridG
    (points,
     bc,
     tess,
#ifdef RICH_MPI
     ParallelHandler(meta),
#endif // RICH_MPI
     NumberIt);
}
