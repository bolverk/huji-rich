#include <cassert>
#include "hdsim_3d.hpp"
#include "../../3D/GeometryCommon/HilbertOrder3D.hpp"
#include "../../misc/utils.hpp"
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
	const SourceTerm3D &source,
	TracerStickerNames &tsn
#ifdef RICH_MPI
	,const ProcessorUpdate3D* proc_update
#endif
) :
	tess_(tess), 
#ifdef RICH_MPI
	tproc_(tproc),
#endif
	eos_(eos), cells_(cells),extensive_(),pm_(pm), tsc_(tsc), fc_(fc), cu_(cu),eu_(eu),source_(source),tsn_(tsn),pt_()
#ifdef RICH_MPI
	,proc_update_(proc_update)
#endif
	, Max_ID_(0)
{
#ifdef RICH_MPI
	int ws = 0, rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	assert(tess.GetPointNo() == cells.size());
	// sort tracers and stickers
	size_t N = tess.GetPointNo();
	vector<size_t> tindex = sort_index(tsn_.tracer_names);
	vector<size_t> sindex = sort_index(tsn_.sticker_names);
	tsn_.tracer_names = VectorValues(tsn_.tracer_names, tindex);
	tsn_.sticker_names = VectorValues(tsn_.sticker_names, sindex);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < tindex.size(); ++j)
			cells_[i].tracers[j] = cells[i].tracers[tindex[j]];
		for (size_t j = 0; j < sindex.size(); ++j)
			cells_[i].stickers[j] = cells[i].stickers[sindex[j]];
	}
	// Is this a new start?
	assert(cells_.size() > 1);
	if (cells_[0].ID == cells_[1].ID)
	{
		size_t nstart = 0;
#ifdef RICH_MPI
		std::vector<size_t> nrecv(static_cast<size_t>(ws), 0);
		size_t nsend = N;
		MPI_Allgather(&nsend, 1, MPI_UNSIGNED_LONG_LONG, &nrecv[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
		for (int i = 0; i < rank; ++i)
			nstart += nrecv[static_cast<size_t>(i)];
#endif
		for (size_t i = 0; i < N; ++i)
			cells_[i].ID = nstart + i;
		Max_ID_ =cells_.back().ID;
#ifdef RICH_MPI
		for (size_t i = static_cast<size_t>(rank + 1); i < static_cast<size_t>(ws); ++i)
			Max_ID_ += nrecv[i];
#endif
	}
	else
	{
		size_t maxid = 0;
		for (size_t i = 0; i < N; ++i)
			maxid = std::max(maxid, cells[i].ID);
#ifdef RICH_MPI
		MPI_Allreduce(&maxid, &Max_ID_, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
		Max_ID_ = maxid;
#endif
	}

#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	extensive_.resize(N);
	for (size_t i = 0; i < N; ++i)
		PrimitiveToConserved(cells_[i], tess.GetVolume(i), extensive_[i]);
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
		{
			res[i] += other[i];
			res[i] *= 0.5;
		}
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
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime(), tsn_);
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel, tsn_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	//dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime(), tsn_);
	vector<Conserved3D> fluxes;
	fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt, tsn_);
	vector<Conserved3D> mid_extensives(extensive_);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, tsn_, mid_extensives);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), tsn_);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D> &mesh = tess_.GetMeshPoints();
		mesh.resize(tess_.GetPointNo());
		vector<size_t> order = HilbertOrder3D(mesh);
		mesh = VectorValues(mesh, order);
		mid_extensives = VectorValues(mid_extensives, order);
		extensive_ = VectorValues(extensive_, order);
		cells_ = VectorValues(cells_, order);
		point_vel = VectorValues(point_vel, order);
	}

#ifdef RICH_MPI
	if (proc_update_ != 0)
		proc_update_->Update(tproc_, tess_);
#endif
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
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, tsn_, mid_extensives);
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
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, tsn_, extensive_);
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

size_t & HDSim3D::GetMaxID(void)
{
	return Max_ID_;
}
