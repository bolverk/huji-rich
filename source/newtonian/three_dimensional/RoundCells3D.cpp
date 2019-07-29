#include "RoundCells3D.hpp"
#include "../../misc/utils.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

RoundCells3D::RoundCells3D(const PointMotion3D& pm, const EquationOfState& eos, Vector3D const& ll, Vector3D const& ur,
	double chi, double eta, bool cold, double min_dw, double dt_speed, vector<std::string> no_move) : pm_(pm), eos_(eos), ll_(ll), ur_(ur), chi_(chi),
	eta_(eta), cold_(cold), min_dw_(min_dw),dt_speed_(dt_speed),no_move_(no_move) {}

namespace
{
	void SlowDown(Vector3D &velocity, Tessellation3D const& tess, double R, size_t index, vector<Vector3D> & velocities,
		vector<char> const& nomove)
	{
		if (nomove[index] == 1)
			return;
		Vector3D const& point = tess.GetMeshPoint(index);
		vector<size_t> neigh = tess.GetNeighbors(index);
		size_t N = neigh.size();
		size_t min_loc = 0;
		double min_d = fastabs(point - tess.GetMeshPoint(neigh[0]));
		for (size_t i = 1; i < N; ++i)
		{
			double d = fastabs(point - tess.GetMeshPoint(neigh[i]));
			if (d < min_d)
			{
				min_d = d;
				min_loc = i;
			}
		}
		if (min_d < 0.3*R)
		{
			size_t face = tess.GetCellFaces(index)[min_loc];
			size_t other = tess.GetFaceNeighbors(face).first == index ? tess.GetFaceNeighbors(face).second :
				tess.GetFaceNeighbors(face).first;
			Vector3D parallel = tess.FaceCM(face) - 0.5*(tess.GetMeshPoint(index) + tess.GetMeshPoint(other));
			parallel *= (1.0 / fastabs(parallel));
			Vector3D v_par = parallel * ScalarProd(velocity, parallel);
			Vector3D v_norm = velocity - parallel * ScalarProd(velocity, parallel);
			if (!tess.IsPointOutsideBox(other))
			{
				if (nomove.at(other) == 1)
				{
					v_par = Vector3D();
				}
				else
				{
					if (index < other)
					{
						v_par *= 0.5;
						Vector3D temp = parallel * ScalarProd(velocities[other], parallel);
						v_par += 0.5*temp;
					}
					else
						v_par = parallel * ScalarProd(velocities[other], parallel);
				}
			}
			velocity = v_norm + v_par;
		}
	}

	void CorrectPointsOverShoot(vector<Vector3D> &v, double dt, Tessellation3D const& tess, Vector3D const& ll,
		Vector3D const& ur)
	{
		// check that we don't go outside grid
		size_t n = tess.GetPointNo();
		const double inv_dt = 1.0 / dt;
		for (size_t i = 0; i < n; ++i)
		{
			Vector3D point(tess.GetMeshPoint(i));
			double R = tess.GetWidth(i);
			if ((v[i].x*dt * 2 + point.x) > ur.x)
			{
				double factor = 0.25*(ur.x - point.x)*inv_dt / fastabs(v[i]);
				if (R*0.1 > (ur.x - point.x))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].y*dt * 2 + point.y) > ur.y)
			{
				double factor = 0.25*(ur.y - point.y)*inv_dt / fastabs(v[i]);
				if (R*0.1 > (ur.y - point.y))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].z*dt * 2 + point.z) > ur.z)
			{
				double factor = 0.25*(ur.z - point.z)*inv_dt / fastabs(v[i]);
				if (R*0.1 > (ur.z - point.z))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].x*dt * 2 + point.x) < ll.x)
			{
				double factor = 0.25*(point.x - ll.x)*inv_dt / fastabs(v[i]);
				if (R*0.1 > (point.x - ll.x))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].y*dt * 2 + point.y) < ll.y)
			{
				double factor = 0.25*(point.y - ll.y)*	inv_dt / fastabs(v[i]);
				if (R*0.1 > (point.y - ll.y))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}
			if ((v[i].z*dt * 2 + point.z) < ll.z)
			{
				double factor = 0.25*(point.z - ll.z)*	inv_dt / fastabs(v[i]);
				if (R*0.1 > (point.z - ll.z))
					factor *= -0.1;
				v[i] = v[i] * factor;
			}

		}
		return;
	}
}

void RoundCells3D::calc_dw(Vector3D &velocity, size_t i, const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	TracerStickerNames const& tracerstickernames, vector<Vector3D> & velocities, vector<char> const& nomove) const
{
	const Vector3D r = tess.GetMeshPoint(i);
	const Vector3D s = tess.GetCellCM(i);
	const double d = fastabs(s - r);
	const double R = tess.GetWidth(i);
	if (d < 0.9*eta_*R)
		return;
	double c = 0;
#ifdef RICH_DEBUG
	try
	{
#endif
		c =std::max(eos_.dp2c(cells[i].density, cells[i].pressure,cells[i].tracers, tracerstickernames.tracer_names), min_dw_);
#ifdef RICH_DEBUG
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Error RoundCells3D::calc_dw", 0);
		eo.AddEntry("Cell number", i);
		throw eo;
	}
#endif
	velocity += chi_ * c*(s - r) / std::max(R, d);
	SlowDown(velocity, tess, R, i, velocities, nomove);
}

void RoundCells3D::calc_dw(Vector3D &velocity, size_t i, const Tessellation3D& tess, double dt, vector<ComputationalCell3D> const& cells,
	TracerStickerNames const& tracerstickernames, vector<Vector3D> & velocities, vector<char> const& nomove)const
{
	const Vector3D r = tess.GetMeshPoint(i);
	const Vector3D s = tess.GetCellCM(i);
	const double d = fastabs(s - r);
	const double R = tess.GetWidth(i);
	if (d < 0.9*eta_*R)
		return;
	vector<size_t> neigh = tess.GetNeighbors(i);
	size_t N = neigh.size();
	double cs = 0;
#ifdef RICH_DEBUG
	try
	{
#endif
		cs = eos_.dp2c(cells[i].density, cells[i].pressure,	cells[i].tracers, tracerstickernames.tracer_names);
#ifdef RICH_DEBUG
	}
	catch (UniversalError &eo)
	{
		eo.AddEntry("Error RoundCells3D::calc_dw", 0);
		eo.AddEntry("Cell number", i);
		throw eo;
	}
#endif
	for (size_t j = 0; j < N; ++j)
	{
		if (tess.IsPointOutsideBox(neigh[j]))
			continue;
#ifdef RICH_DEBUG
		try
		{
#endif
			cs = std::max(cs, eos_.dp2c(cells[neigh[j]].density, cells[neigh[j]].pressure,
				cells[static_cast<size_t>(neigh[j])].tracers, tracerstickernames.tracer_names));
			cs = std::max(cs, fastabs(cells[neigh[j]].velocity-cells[i].velocity));
#ifdef RICH_DEBUG
			if (!std::isfinite(cs))
				throw UniversalError("Bad cs in roundcells");
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error RoundCells3D::calc_dw", 0);
			eo.AddEntry("Neigh Cell number", neigh[j]);
			eo.AddEntry("Neigh index", j);
			throw eo;
	}
#endif
		
}
	const double c_dt = std::max(std::max(dt_speed_*d / dt, cs), min_dw_);
	velocity += chi_ * c_dt*(s - r) / std::max(R, d);
	SlowDown(velocity, tess, R, i, velocities, nomove);
}

void RoundCells3D::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
	double time, TracerStickerNames const& tracerstickernames, vector<Vector3D> &res) const
{
	pm_(tess, cells, time, tracerstickernames, res);
	const size_t n = tess.GetPointNo();
	if (n == 0)
		return;
	size_t Nstick = tracerstickernames.sticker_names.size();
	vector<char> nomove(n, 0);
	vector<size_t> no_move_indeces;
	for (size_t i = 0; i < Nstick; ++i)
	{
		vector<std::string>::const_iterator it = std::find(no_move_.begin(), no_move_.end(),
			tracerstickernames.sticker_names.at(i));
		if (it != no_move_.end())
			no_move_indeces.push_back(i);
	}
	Nstick = no_move_indeces.size();
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < Nstick; ++j)
		{
			if (cells[i].stickers[no_move_indeces[j]])
			{
				res[i] = Vector3D();
				nomove[i] = 1;
				break;
			}
		}
	}
#ifndef RICH_MPI
	for (size_t i = 0; i < n; ++i)
		SlowDown(res[i], tess, tess.GetWidth(i), i, res, nomove);
#endif
}

void RoundCells3D::ApplyFix(Tessellation3D const& tess, vector<ComputationalCell3D> const& cells, double time,
	double dt, vector<Vector3D> &velocities, TracerStickerNames const& tracerstickernames)const
{
	pm_.ApplyFix(tess, cells, time, dt, velocities, tracerstickernames);
#ifdef RICH_MPI
	Vector3D vdummy;
	MPI_exchange_data(tess, velocities, true,&vdummy);
#endif
	const size_t n = tess.GetPointNo();
	/*if (n == 0)
		return;*/
	vector<char> nomove(n, 0);
	size_t Nstick = tracerstickernames.sticker_names.size();
	vector<size_t> no_move_indeces;
	for (size_t i = 0; i < Nstick; ++i)
	{
		vector<std::string>::const_iterator it = std::find(no_move_.begin(), no_move_.end(),
			tracerstickernames.sticker_names.at(i));
		if (it != no_move_.end())
			no_move_indeces.push_back(i);
	}
	Nstick = no_move_indeces.size();
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < Nstick; ++j)
		{
			if (cells[i].stickers[no_move_indeces[j]])
			{
				velocities[i] = Vector3D();
				nomove[i] = 1;
				break;
			}
		}
	}
#ifdef RICH_MPI
	MPI_exchange_data(tess, nomove, true);
#endif

	if (cold_)
	{
		for (size_t i = 0; i < n; ++i)
		{
			if (nomove[i] == 0)
				calc_dw(velocities.at(i), i, tess, dt, cells, tracerstickernames, velocities, nomove);
			else
				velocities.at(i) = Vector3D();
		}
	}
	else
	{
		for (size_t i = 0; i < n; ++i)
		{
			if (nomove[i] == 0)
				calc_dw(velocities.at(i), i, tess, cells, tracerstickernames, velocities, nomove);
			else
				velocities.at(i) = Vector3D();
		}
	}
#ifdef RICH_MPI
	for (size_t i = 0; i < n; ++i)
		SlowDown(velocities[i], tess, tess.GetWidth(i), i, velocities, nomove);
#endif
	velocities.resize(n);
	CorrectPointsOverShoot(velocities, dt, tess, ll_, ur_);
	}

void RoundCells3D::ChangeBox(Vector3D const & ll, Vector3D const & ur)
{
	ll_ = ll;
	ur_ = ur;
}
