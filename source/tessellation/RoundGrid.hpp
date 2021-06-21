/*! \file RoundGrid.hpp
  \brief Makes the initial cells rounder
  \author Elad Steinberg
 */

#ifndef ROUNDGRID
#define ROUNDGRID 1
#include "VoronoiMesh.hpp"

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

  virtual ~TessellationHandler(void);
};

class SerialHandler: public TessellationHandler
{
public:

  void initialise
  (Tessellation& tess,
   const vector<Vector2D>& points,
   const OuterBoundary& bc) const;

  void update
  (Tessellation& tess, 
   const vector<Vector2D>& points) const;
};

#ifdef RICH_MPI
class ParallelHandler: public TessellationHandler
{
public:

  ParallelHandler(const Tessellation& meta);

  void initialise
  (Tessellation& tess,
   const vector<Vector2D>& points,
   const OuterBoundary& bc) const;

  void update
  (Tessellation& tess,
   const vector<Vector2D>& points) const;

private:
  const Tessellation& meta_;
};
#endif // RICH_MPI

vector<Vector2D> RoundGrid
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
#ifdef RICH_MPI
 const TessellationHandler& th,
#endif
 int NumberIt=10);

/*!
	\brief Makes the cells rounder
	\param points The initial points
	\param bc The outer boundary conditions
	\param NumberIt The number of correction iterations
	\param tess The tessellation
*/
#ifdef RICH_MPI
//!	\param tproc The tessellation of processors
#endif
/*!
	\return The points that give a rounder tessellation
*/
vector<Vector2D> RoundGrid
(vector<Vector2D> const& points,
 const OuterBoundary& bc,
 Tessellation& tess,
#ifdef RICH_MPI
 const TessellationHandler& th,
#endif
 int NumberIt=10);

#endif //ROUNDGRID
