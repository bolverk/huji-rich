#include <cassert>
#include "hdsim_3d.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

Tessellation3D& HDSim3D::getTesselation(void)
{
	return tess_;
}

#ifdef RICH_MPI
const Tessellation3D& HDSim3D::getProcTesselation(void) const
{
	return tproc_;
}
#endif

vector<ComputationalCell3D>& HDSim3D::getCells(void)
{
	return cells_;
}

vector<Conserved3D>& HDSim3D::getExtensives(void)
{
	return extensive_;
}


const vector<Conserved3D>& HDSim3D::getExtensives(void) const
{
	return extensive_;
}

HDSim3D::ProgressTracker::ProgressTracker(void) :
	time(0), cycle(0) {}

void HDSim3D::ProgressTracker::update(double dt)
{
	++cycle;
	time += dt;
}

double HDSim3D::ProgressTracker::getTime(void) const
{
	return time;
}

size_t HDSim3D::ProgressTracker::getCycle(void) const
{
	return cycle;
}

HDSim3D::HDSim3D(Tessellation3D& tess,
#ifdef RICH_MPI
	Tessellation3D& tproc,
#endif//RICH_MPI
	const vector<ComputationalCell3D>& cells,
	const EquationOfState& eos,
	const PointMotion3D& pm,
	const TimeStepFunction3D& tsc,
	const FluxCalculator3D& fc,
	const CellUpdater3D& cu,
	const ExtensiveUpdater3D& eu,
	const TracerStickerNames tsn) :
	tess_(tess), 
#ifdef RICH_MPI
	tproc_(tproc),
#endif
	eos_(eos), cells_(cells),extensive_(),pm_(pm), tsc_(tsc), fc_(fc), cu_(cu),eu_(eu),tsn_(tsn),pt_()
{
	assert(tess.GetPointNo() == cells.size());
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	size_t N = tess.GetPointNo();
	extensive_.resize(N);
	for (size_t i = 0; i < N; ++i)
		PrimitiveToConserved(cells_[i], tess.GetVolume(i), extensive_[i], eos_, tsn_);
}

namespace
{
	void CalcFaceVelocities(Tessellation3D const& tess, vector<Vector3D> const& point_vel, vector<Vector3D> &res)
	{
		size_t N = tess.GetTotalFacesNumber();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
			if (tess.BoundaryFace(i))
				res[i] = Vector3D();
			else
				res[i]=tess.CalcFaceVelocity(i, point_vel[tess.GetFaceNeighbors(i).first], point_vel[tess.GetFaceNeighbors(i).second]);
	}

	void UpdateTessellation(Tessellation3D &tess, vector<Vector3D> &point_vel,double dt
#ifdef RICH_MPI
		,Tessellation3D const& tproc
#endif
		)
	{
		vector<Vector3D> points = tess.GetMeshPoints();
		size_t N = tess.GetPointNo();
		points.resize(N);
		for (size_t i = 0; i < N; ++i)
			points[i] += point_vel[i] * dt;
		tess.Build(points
#ifdef RICH_MPI
			,tproc
#endif
			);
	}

	void ExtensiveAvg(vector<Conserved3D> &res, vector<Conserved3D> const& other)
	{
		assert(res.size() == other.size());
		size_t N = res.size();
		for (size_t i = 0; i < N; ++i)
			res[i] = 0.5*(res[i] + other[i]);
	}
}


void HDSim3D::timeAdvance2(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), tsn_, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	const double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime(), tsn_);
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	vector<Conserved3D> fluxes;
	fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt, tsn_);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), tsn_);
	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		,tproc_
#endif
		);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif

	cu_(cells_, eos_, tess_, mid_extensives, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif

	pt_.update(dt);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt, tsn_);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), tsn_);
	ExtensiveAvg(extensive_, mid_extensives);
	cu_(cells_, eos_, tess_, extensive_, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
}

void HDSim3D::timeAdvance(void)
{
	vector<Vector3D> point_vel,face_vel;
	pm_(tess_, cells_, pt_.getTime(), tsn_, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	const double dt = tsc_(tess_, cells_, eos_,face_vel,pt_.getTime(),tsn_);
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	vector<Conserved3D> fluxes;
	fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt, tsn_);
	eu_(fluxes, tess_, dt, cells_, extensive_, pt_.getTime(), tsn_);
	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
		);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
#endif
	cu_(cells_, eos_, tess_, extensive_, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	pt_.update(dt);
}

const Tessellation3D& HDSim3D::getTesselation(void) const
{
	return tess_;
}

const vector<ComputationalCell3D>& HDSim3D::getCells(void) const
{
	return cells_;
}

double HDSim3D::GetTime(void)const
{
	return pt_.getTime();
}

TracerStickerNames HDSim3D::GetTracerStickerNames(void)const
{
	return tsn_;
}

size_t HDSim3D::GetCycle(void)const
{
	return static_cast<size_t>(pt_.getCycle());
}

void HDSim3D::SetCycle(size_t cycle)
{
	pt_.cycle = cycle;
}

void HDSim3D::SetTime(double t)
{
	pt_.time = t;
}
