#include "hold_still.hpp"

HoldStill::Condition::~Condition(void) {}

HoldStill::HoldStill(PointMotion& raw,
		     const Condition& cond):
  raw_(raw), cond_(cond) {}

Vector2D HoldStill::CalcVelocity(int index,
				 Tessellation const& tess,
				 vector<Primitive> const& cells,
				 double time)
{
  if(cond_(index,tess,cells,time))
    return Vector2D(0,0);
  else
    return raw_.CalcVelocity(index,tess,cells,time);
}

vector<Vector2D> HoldStill::calcAllVelocities(Tessellation const& tess,
					      vector<Primitive> const& cells,
					      double time,
					      vector<CustomEvolution*> &cevolve,
					      const vector<vector<double> >& tracers)
{
  vector<Vector2D> res = raw_.calcAllVelocities(tess,cells,time,cevolve,tracers);
  for(size_t i=0;i<res.size();++i)
  {
    if(cond_(static_cast<int>(i),tess,cells,time))
      res[i] = Vector2D(0,0);
  }
  return res;
}

void HoldStill::ApplyFix(Tessellation const& tess, vector<Primitive> const& cells, double time,
	vector<CustomEvolution*> &cevolve, const vector<vector<double> >& tracers, double dt, vector < Vector2D >
	& velocities)
{
	raw_.ApplyFix(tess, cells, time, cevolve, tracers, dt, velocities);
	const size_t n = static_cast<size_t>(tess.GetPointNo());
	for (size_t i = 0; i<n; ++i)
	{
		if (cond_(static_cast<int>(i), tess, cells, time))
			velocities[i] = Vector2D(0, 0);
	}
}