#include <cassert>
#include "hdsim_3d.hpp"
#include "extensive_generator.hpp"

HDSim3D::ProgressTracker::ProgressTracker(void):
  time_(0), cycle_(0) {}

void HDSim3D::ProgressTracker::update(double dt)
{
  ++cycle_;
  time_ += dt;
}

double HDSim3D::ProgressTracker::getTime(void) const
{
  return time_;
}

double HDSim3D::ProgressTracker::getCycle(void) const
{
  return cycle_;
}

HDSim3D::HDSim3D(const Tessellation3D& tess,
		 const vector<ComputationalCell>& cells,
		 const EquationOfState& eos,
		 const PointMotion3D& pm,
		 const TimeStepCalculator& tsc,
		 const FluxCalculator& fc,
		 const CellUpdater& cu):
  tess_(tess), eos_(eos), cells_(cells),
  extensive_(serial_generate(ExtensiveGenerator(cells,tess,eos))),
  pm_(pm), tsc_(tsc), fc_(fc), cu_(cu)
{
  assert(tess.GetPointNo()==cells.size());
}

namespace {
  class PointVelocitiesCalculator: public Index2Member<Vector3D>
  {
  public:

    PointVelocitiesCalculator(const PointMotion3D& pm,
			      const Tessellation3D& tess):
      pm_(pm), tess_(tess) {}

    size_t getLength(void) const
    {
      return tess_.GetPointNo();
    }

    Vector3D operator()(size_t i) const
    {
      return pm_(tess_.GetMeshPoint(i));
    }

  private:
    const PointMotion3D& pm_;
    const Tessellation3D& tess_;
  };

  void update_extensive(const vector<Conserved3D>& fluxes,
			double dt,
			const Tessellation3D& tess,
			vector<Conserved3D>& extensive)
  {
    for(size_t i=0;i<tess.GetTotalFacesNumber();++i){
      const Conserved3D delta = dt*tess.GetEdge(i).GetArea()*fluxes[i];
      if(tess.GetEdge(i).neighbors.first.isInDomain())
	extensive[tess.GetEdge(i).neighbors.first.getIndex()] -= delta;
      if(tess.GetEdge(i).neighbors.second.isInDomain())
	extensive[tess.GetEdge(i).neighbors.second.getIndex()] += delta;
    }
  }

  class AllCellsUpdater: public Index2Member<ComputationalCell>
  {
  public:

    AllCellsUpdater(const Tessellation3D& tess,
		    const vector<Conserved3D>& extensive,
		    const EquationOfState& eos,
		    const CellUpdater& cu):
      tess_(tess), extensive_(extensive), eos_(eos), cu_(cu) {}

    size_t getLength(void) const
    {
      return extensive_.size();
    }

    ComputationalCell operator()(size_t i) const
    {
      return cu_(extensive_[i]/tess_.GetVolume(i),eos_);
    }

  private:
    const Tessellation3D& tess_;
    const vector<Conserved3D>& extensive_;
    const EquationOfState& eos_;
    const CellUpdater& cu_;
  };
}

void HDSim3D::timeAdvance(void)
{
  const double dt = tsc_(tess_,cells_,eos_);
  const vector<Vector3D> point_velocities = 
    serial_generate(PointVelocitiesCalculator(pm_,tess_));
  update_extensive(fc_(tess_,cells_,point_velocities),
		   dt,tess_,extensive_);
  cells_ = serial_generate(AllCellsUpdater(tess_,extensive_,eos_,cu_));  
  
}
