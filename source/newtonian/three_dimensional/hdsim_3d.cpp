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

Tessellation3D& HDSim3D::getProcTesselation(void)
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

void HDSim3D::ProgressTracker::updateTime(double dt)
{
	time += dt;
}

void HDSim3D::ProgressTracker::updateCycle()
{
	++cycle;
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
	const SourceTerm3D& source,
		 const pair<vector<string>, vector<string> >& tsn,
	bool SR
#ifdef RICH_MPI
	, const ProcessorUpdate3D* proc_update
#endif
	, bool new_start
#ifdef RICH_MPI
		 , const double maxload
#endif // RICH_MPI
) :
	tess_(tess),
#ifdef RICH_MPI
	tproc_(tproc),
#endif
	eos_(eos), cells_(cells), extensive_(), pm_(pm), tsc_(tsc), fc_(fc), cu_(cu), eu_(eu), source_(source), pt_()
#ifdef RICH_MPI
	, proc_update_(proc_update)
#endif
	, Max_ID_(0)
#ifdef RICH_MPI
	, maxload_(maxload)
#endif // RICH_MPI
{
#ifdef RICH_MPI
	int ws = 0, rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	const bool validity_check = tess.GetPointNo() == cells.size();
	assert(validity_check);
	assert(tsn.second.size() <= MAX_STICKERS);
	assert(tsn.first.size() <= MAX_TRACERS);
	// sort tracers and stickers
	size_t N = tess.GetPointNo();
	vector<size_t> tindex = sort_index(tsn.first);
	vector<size_t> sindex = sort_index(tsn.second);
	ComputationalCell3D::tracerNames = VectorValues(tsn.first, tindex);
	ComputationalCell3D::stickerNames = VectorValues(tsn.second, sindex);
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < tindex.size(); ++j)
			cells_[i].tracers[j] = cells[i].tracers[tindex[j]];
		for (size_t j = 0; j < sindex.size(); ++j)
			cells_[i].stickers[j] = cells[i].stickers[sindex[j]];
	}
	// Is this a new start?
	if (new_start)
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
		Max_ID_ = nstart + N - 1;
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
	ComputationalCell3D cdummy;
	MPI_exchange_data(tess_, cells_, true, &cdummy);
#endif
	extensive_.resize(N);
	if (SR)
	{
		for (size_t i = 0; i < N; ++i)
			PrimitiveToConservedSR(cells_[i], tess.GetVolume(i), extensive_[i], eos_);
	}
	else
	{
		for (size_t i = 0; i < N; ++i)
			PrimitiveToConserved(cells_[i], tess.GetVolume(i), extensive_[i]);
	}
}

namespace
{
	void CalcFaceVelocities(Tessellation3D const& tess, vector<Vector3D> const& point_vel, vector<Vector3D>& res)
	{
		size_t N = tess.GetTotalFacesNumber();
		res.resize(N);
		for (size_t i = 0; i < N; ++i)
		{
			if (tess.BoundaryFace(i))
				res[i] = Vector3D();
			else
			{
				try
				{
					res[i] = tess.CalcFaceVelocity(i, point_vel[tess.GetFaceNeighbors(i).first], point_vel[tess.GetFaceNeighbors(i).second]);
				}
				catch (UniversalError & /*eo*/)
				{
					throw;
				}
			}
		}
	}

	void UpdateTessellation(Tessellation3D& tess, const vector<Vector3D>& point_vel, double dt
#ifdef RICH_MPI
		, Tessellation3D const& tproc
#endif
		, std::vector<Vector3D> const* orgpoints = nullptr)
	{
		vector<Vector3D> points;
		if (orgpoints == nullptr)
			points = tess.getMeshPoints();
		else
			points = *orgpoints;
		points.resize(tess.GetPointNo());
		size_t N = points.size();
		for (size_t i = 0; i < N; ++i)
			points[i] += point_vel[i] * dt;
		tess.Build(points
#ifdef RICH_MPI
			, tproc
#endif
		);
	}

	void ExtensiveAvg(vector<Conserved3D>& res, vector<Conserved3D> const& other)
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
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	Vector3D vdummy;
	MPI_exchange_data(tess_, point_vel, true, &vdummy);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true, &vdummy);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values = 
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D>& mesh = tess_.accessMeshPoints();
		mesh.resize(tess_.GetPointNo());
		vector<size_t> order = HilbertOrder3D(mesh);
		mesh = VectorValues(mesh, order);
		mid_extensives = VectorValues(mid_extensives, order);
		extensive_ = VectorValues(extensive_, order);
		cells_ = VectorValues(cells_, order);
		point_vel = VectorValues(point_vel, order);
	}

#ifdef RICH_MPI
	int ntotal = 0;
	double load = 6.0;
	ComputationalCell3D cdummy;
	Conserved3D edummy;
	while (load > maxload_)
	{
		if (proc_update_ != 0)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			int rank = 0;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			vector<size_t> selfindex;
			vector<vector<size_t> > sentpoints;
			vector<int> sentproc;
			MPI_Barrier(MPI_COMM_WORLD);
			proc_update_->Update(tproc_, tess_);
			std::vector<Vector3D> &oldpoints = tess_.accessMeshPoints();
			oldpoints.resize(tess_.GetPointNo());
			std::vector<Vector3D> newpoints = tess_.UpdateMPIPoints(tproc_, rank, oldpoints, selfindex, sentproc, sentpoints);
			MPI_Barrier(MPI_COMM_WORLD);
			tess_.GetPointNo() = newpoints.size();
			tess_.GetSentPoints() = sentpoints;
			tess_.GetSentProcs() = sentproc;
			tess_.GetSelfIndex() = selfindex;
			tess_.accessMeshPoints() = newpoints;
			// Keep relevant points
			MPI_exchange_data(tess_, mid_extensives, false, &edummy);
			MPI_exchange_data(tess_, extensive_, false, &edummy);
			MPI_exchange_data(tess_, cells_, false, &cdummy);
			MPI_exchange_data(tess_, point_vel, false, &vdummy);
			load = proc_update_->GetLoadImbalance(tess_, ntotal);
		}
		else
			load = 0.0;
	}
#endif
	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
	);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false, &edummy);
	MPI_exchange_data(tess_, extensive_, false, &edummy);
	MPI_exchange_data(tess_, cells_, false, &cdummy);
	MPI_exchange_data(tess_, point_vel, false, &vdummy);
	MPI_exchange_data(tess_, point_vel, true, &vdummy);
#endif

cu_(cells_, eos_, tess_, mid_extensives);
#ifdef RICH_MPI
MPI_exchange_data(tess_, cells_, true, &cdummy);
#endif

pt_.updateTime(dt);
pt_.updateCycle();
CalcFaceVelocities(tess_, point_vel, face_vel);
face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
ExtensiveAvg(extensive_, mid_extensives);
cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
MPI_exchange_data(tess_, cells_, true, &cdummy);
#endif
}

void HDSim3D::timeAdvance(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	Vector3D vdummy;
	MPI_exchange_data(tess_, point_vel, true, &vdummy);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	const double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true, &vdummy);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values = 
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, extensive_);
	eu_(fluxes, tess_, dt, cells_, extensive_, pt_.getTime(), face_vel, face_values);
	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
	);
#ifdef RICH_MPI
	// Keep relevant points
	ComputationalCell3D cdummy;
	Conserved3D edummy;
	MPI_exchange_data(tess_, extensive_, false, &edummy);
	MPI_exchange_data(tess_, cells_, false, &cdummy);
#endif
	cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true, &cdummy);
#endif
	pt_.updateTime(dt);
	pt_.updateCycle();
}


void HDSim3D::timeAdvance3(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values = 
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), 0.5 * dt);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, 0.5 * dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), 0.5 * dt, mid_extensives);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D>& mesh = tess_.accessMeshPoints();
		mesh.resize(tess_.GetPointNo());
		vector<size_t> order = HilbertOrder3D(mesh);
		mesh = VectorValues(mesh, order);
		mid_extensives = VectorValues(mid_extensives, order);
		extensive_ = VectorValues(extensive_, order);
		cells_ = VectorValues(cells_, order);
		point_vel = VectorValues(point_vel, order);
		//du1 = VectorValues(du1, order);
	}
	std::vector<Vector3D> oldpoints = tess_.accessMeshPoints();
	oldpoints.resize(tess_.GetPointNo());
#ifdef RICH_MPI
	if (proc_update_ != 0)
		proc_update_->Update(tproc_, tess_);
#endif
	UpdateTessellation(tess_, point_vel, 0.5 * dt
#ifdef RICH_MPI
		, tproc_
#endif
	);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	//MPI_exchange_data(tess_, du1, false);
	MPI_exchange_data(tess_, oldpoints, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	std::vector<Conserved3D> u1 = mid_extensives;
	cu_(cells_, eos_, tess_, mid_extensives);


#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif

	pt_.updateTime(0.5 * dt);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), 2 * dt);
	//mid_extensives = extensive_;
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), 2 * dt,  mid_extensives);
	eu_(fluxes, tess_, 2 * dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	mid_extensives = mid_extensives - 3 * (u1 - extensive_);

	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
		, &oldpoints);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, u1, false);
	//MPI_exchange_data(tess_, du2, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	std::vector<Conserved3D> u2 = mid_extensives;
	cu_(cells_, eos_, tess_, mid_extensives);


#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	pt_.updateTime(0.5 * dt);
	pt_.updateCycle();
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt / 6);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt / 6,  mid_extensives);
	eu_(fluxes, tess_, dt / 6, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	extensive_ = mid_extensives - (1.0 / 3.0) * (2 * u2 + extensive_) + u1;
	cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
}

void HDSim3D::timeAdvance33(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values = 
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), dt);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D>& mesh = tess_.accessMeshPoints();
		mesh.resize(tess_.GetPointNo());
		vector<size_t> order = HilbertOrder3D(mesh);
		mesh = VectorValues(mesh, order);
		mid_extensives = VectorValues(mid_extensives, order);
		extensive_ = VectorValues(extensive_, order);
		cells_ = VectorValues(cells_, order);
		point_vel = VectorValues(point_vel, order);
	}
	std::vector<Vector3D> oldpoints = tess_.accessMeshPoints();
	oldpoints.resize(tess_.GetPointNo());
#ifdef RICH_MPI
	if (proc_update_ != 0)
		proc_update_->Update(tproc_, tess_);
#endif
	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
	);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	MPI_exchange_data(tess_, oldpoints, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	cu_(cells_, eos_, tess_, mid_extensives);


#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif

	pt_.updateTime(dt);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
	//mid_extensives = extensive_;
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	mid_extensives = 0.25 * mid_extensives + 0.75 * extensive_;

	UpdateTessellation(tess_, point_vel, dt / 2
#ifdef RICH_MPI
		, tproc_
#endif
		, &oldpoints);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, oldpoints, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	cu_(cells_, eos_, tess_, mid_extensives);


#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	pt_.updateTime(-0.5 * dt);
	pt_.updateCycle();
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	extensive_ = 0.33333333333333333333333 * (2 * mid_extensives + extensive_);

	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
		, &oldpoints);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
#endif

	cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
}

void HDSim3D::timeAdvance32(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values =
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), 0.5 * dt);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D>& mesh = tess_.accessMeshPoints();
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
		, tproc_
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
	std::vector<Conserved3D> u1 = mid_extensives;
	cu_(cells_, eos_, tess_, mid_extensives);

#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif

	pt_.updateTime(dt);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	mid_extensives = 0.5 * (mid_extensives + extensive_);
	cu_(cells_, eos_, tess_, mid_extensives);

#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	pt_.updateCycle();
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	//extensive_ = 0.333333333333333333*(extensive_ + u1 + mid_extensives);
	extensive_ = 0.333333333333333333 * (extensive_ + u1 + mid_extensives);
	cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
}

void HDSim3D::timeAdvance4(void)
{
	vector<Vector3D> point_vel, face_vel;
	pm_(tess_, cells_, pt_.getTime(), point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif

	CalcFaceVelocities(tess_, point_vel, face_vel);
	double dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	pm_.ApplyFix(tess_, cells_, pt_.getTime(), dt, point_vel);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, point_vel, true);
#endif
	CalcFaceVelocities(tess_, point_vel, face_vel);
	dt = tsc_(tess_, cells_, eos_, face_vel, pt_.getTime());
	vector<Conserved3D> fluxes;
	std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > face_values =
		fc_(fluxes, tess_, face_vel, cells_, extensive_, eos_, pt_.getTime(), 0.5 * dt);
	vector<Conserved3D> mid_extensives(extensive_);
	eu_(fluxes, tess_, 0.5 * dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), 0.5 * dt, mid_extensives);

	if (pt_.cycle % 10 == 0)
	{
		vector<Vector3D>& mesh = tess_.accessMeshPoints();
		mesh.resize(tess_.GetPointNo());
		vector<size_t> order = HilbertOrder3D(mesh);
		mesh = VectorValues(mesh, order);
		mid_extensives = VectorValues(mid_extensives, order);
		extensive_ = VectorValues(extensive_, order);
		cells_ = VectorValues(cells_, order);
		point_vel = VectorValues(point_vel, order);
		//du1 = VectorValues(du1, order);
	}
	std::vector<Vector3D> oldpoints = tess_.accessMeshPoints();
	oldpoints.resize(tess_.GetPointNo());
#ifdef RICH_MPI
	if (proc_update_ != 0)
		proc_update_->Update(tproc_, tess_);
#endif
	UpdateTessellation(tess_, point_vel, 0.5 * dt
#ifdef RICH_MPI
		, tproc_
#endif
	);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	//MPI_exchange_data(tess_, du1, false);
	MPI_exchange_data(tess_, oldpoints, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	cu_(cells_, eos_, tess_, mid_extensives);
	std::vector<Conserved3D> du1 = mid_extensives - extensive_;

#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif

	pt_.updateTime(0.5 * dt);
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), 0.5 * dt);
	//mid_extensives = extensive_;
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), 0.5 * dt, mid_extensives);
	mid_extensives = mid_extensives - du1;
	eu_(fluxes, tess_, 0.5 * dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	cu_(cells_, eos_, tess_, mid_extensives);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	std::vector<Conserved3D> du2 = mid_extensives - extensive_;

	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt, mid_extensives);
	mid_extensives = mid_extensives - du2;
	eu_(fluxes, tess_, dt, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);

	UpdateTessellation(tess_, point_vel, dt
#ifdef RICH_MPI
		, tproc_
#endif
		, &oldpoints);
#ifdef RICH_MPI
	// Keep relevant points
	MPI_exchange_data(tess_, mid_extensives, false);
	MPI_exchange_data(tess_, du1, false);
	MPI_exchange_data(tess_, du2, false);
	//MPI_exchange_data(tess_, du3, false);
	MPI_exchange_data(tess_, extensive_, false);
	MPI_exchange_data(tess_, cells_, false);
	MPI_exchange_data(tess_, point_vel, false);
	MPI_exchange_data(tess_, point_vel, true);
#endif
	cu_(cells_, eos_, tess_, mid_extensives);
	std::vector<Conserved3D> du3 = mid_extensives - extensive_;

#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
	pt_.updateTime(0.5 * dt);
	pt_.updateCycle();
	CalcFaceVelocities(tess_, point_vel, face_vel);
	face_values = fc_(fluxes, tess_, face_vel, cells_, mid_extensives, eos_, pt_.getTime(), dt / 6);
	source_(tess_, cells_, fluxes, point_vel, pt_.getTime(), dt / 6,  mid_extensives);
	mid_extensives = mid_extensives - du3;
	eu_(fluxes, tess_, dt / 6, cells_, mid_extensives, pt_.getTime(), face_vel, face_values);
	extensive_ = mid_extensives + (1.0 / 6.0) * (2 * du1 + 4 * du2 + 2 * du3);
	cu_(cells_, eos_, tess_, extensive_);
#ifdef RICH_MPI
	MPI_exchange_data(tess_, cells_, true);
#endif
}

const Tessellation3D& HDSim3D::getTesselation(void) const
{
	return tess_;
}

const vector<ComputationalCell3D>& HDSim3D::getCells(void) const
{
	return cells_;
}

double HDSim3D::getTime(void)const
{
	return pt_.getTime();
}

size_t HDSim3D::getCycle(void)const
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

size_t& HDSim3D::GetMaxID(void)
{
	return Max_ID_;
}

double HDSim3D::RadiationTimeStep(double const dt, CG::MatrixBuilder const& matrix_builder, std::string const& key)
{
	int total_iters = 0;
	double const CG_eps = 1e-10;
	std::vector<double> new_Er = CG::conj_grad_solver(CG_eps, total_iters, tess_, cells_ , key, dt, matrix_builder);
	size_t const N = tess_.GetPointNo();
	double const max_Er = *std::max_element(new_Er.begin(), new_Er.end());
	double max_diff = 0;
	for(size_t i = 0; i < N; ++i)
	{
		double const diff = std::abs(new_Er[i] - cells_[i].tracers[0] * cells_[i].density) / (new_Er[i] + 0.001 * max_Er);
		max_diff = std::max(max_diff, diff);
	}
	matrix_builder.PostCG(tess_, extensive_, dt, cells_, new_Er, key);
#ifdef RICH_MPI
	ComputationalCell3D cdummy;
	MPI_exchange_data(tess_, cells_, true, &cdummy);
	MPI_Allreduce(MPI_IN_PLACE, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
	pt_.updateTime(dt);
	pt_.updateCycle();
	return dt * std::min(0.01 / max_diff, 1.1);
}