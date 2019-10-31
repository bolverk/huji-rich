#include "ConditionExtensiveUpdater3D.hpp"
#include "../../misc/utils.hpp"
#include <iostream>
#include <cfloat>
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

ConditionExtensiveUpdater3D::Condition3D::~Condition3D() {}

ConditionExtensiveUpdater3D::Action3D::~Action3D() {}

ConditionExtensiveUpdater3D::~ConditionExtensiveUpdater3D() {}

ConditionExtensiveUpdater3D::ConditionExtensiveUpdater3D(const vector<pair<const Condition3D*, const Action3D*> >& sequence) :
	sequence_(sequence) {}

void ConditionExtensiveUpdater3D::operator()(const vector<Conserved3D>& fluxes, const Tessellation3D& tess,
	const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double time,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = tess.GetPointNo();
	std::vector<double> oldEk(N, 0), oldEtherm(N, 0), oldE(N, 0);
	for (size_t i = 0; i < N; ++i)
	{
		oldEk[i] = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		oldEtherm[i] = extensives[i].internal_energy;
		oldE[i] = extensives[i].energy;
	}
	bool entropy = !(std::find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(), std::string("Entropy")) ==
		tracerstickernames.tracer_names.end());
	size_t entropy_index = static_cast<size_t>(std::find(tracerstickernames.tracer_names.begin(),
		tracerstickernames.tracer_names.end(), std::string("Entropy")) - tracerstickernames.tracer_names.begin());
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
	for (size_t i = 0; i < Nfluxes; ++i)
	{
		delta = fluxes[i] * dt * tess.GetArea(i);
		delta.internal_energy = 0;
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;
		if (n0 < N)
		{
			extensives[n0] -= delta;
			extensives[n0].internal_energy -= delta.energy - ScalarProd(cells[n0].velocity, delta.momentum) +
				0.5 * ScalarProd(cells[n0].velocity, cells[n0].velocity) * delta.mass;
		}
		if (n1 < N)
		{
			extensives[n1] += delta;
			extensives[n1].internal_energy += delta.energy - ScalarProd(cells[n1].velocity, delta.momentum) +
				0.5 * ScalarProd(cells[n1].velocity, cells[n1].velocity) * delta.mass;
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		double dEtherm = extensives[i].internal_energy - oldEtherm[i];
		double Eknew = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		double dEk = Eknew - oldEk[i];
		double dE = extensives[i].energy - oldE[i];
		if (dEtherm * (dE - dEk) > 0)
		{
			if (std::abs(dEtherm) > 0.95 * std::abs(dE - dEk) && std::abs(dEtherm) < 1.05 * std::abs(dE - dEk))
				extensives[i].internal_energy = extensives[i].energy - Eknew;
			else
				extensives[i].energy = extensives[i].internal_energy + Eknew;
		}
		else
			extensives[i].energy = extensives[i].internal_energy + Eknew;
		for (size_t j = 0; j < sequence_.size(); ++j)
		{
			if (sequence_[j].first->operator()(i, tess, cells, time, tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, tess, dt, cells, extensives, i, time, tracerstickernames);
				break;
			}
		}
		// check cell
		if (!(extensives[i].mass > 0) || !(extensives[i].energy > 0) || (!(extensives[i].internal_energy > 0) && (!entropy)) ||
			(!std::isfinite(fastabs(extensives[i].momentum))) || (entropy && extensives[i].tracers[entropy_index] < 0))
		{
			int rank = 0;
#ifdef RICH_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
			std::cout << "Bad cell in ExtensiveUpdater3D, cell " << i << " rank " << rank << std::endl;
			std::cout << "mass " << extensives[i].mass << " energy " << extensives[i].energy << " internalE " <<
				extensives[i].internal_energy << " momentum" << abs(extensives[i].momentum) << " volume " << tess.GetVolume(i)
				<< std::endl;
			std::cout << "Old cell, density " << cells[i].density << " pressure " << cells[i].pressure << " vx " <<
				cells[i].velocity.x << " vy " << cells[i].velocity.y << " vz " << cells[i].velocity.z << std::endl;
			for (size_t j = 0; j < tracerstickernames.tracer_names.size(); ++j)
			{
				std::cout << tracerstickernames.tracer_names[j] << " old cell " << cells[i].tracers[j] << 
					" extensive "<<extensives[i].tracers[j]<<std::endl;
			}
			face_vec temp = tess.GetCellFaces(i);
			Conserved3D old_ext(extensives[i]);
			// recover old_extensive
			for (size_t j = 0; j < temp.size(); ++j)
			{
				double Area = tess.GetArea(temp[j]) * dt;
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				if (N1 == i)
				{
					double newEk = 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
					old_ext.mass -= Area * fluxes[temp[j]].mass;
					old_ext.momentum -= Area * fluxes[temp[j]].momentum;
					old_ext.energy -= Area * fluxes[temp[j]].energy;
					old_ext.internal_energy -= Area * fluxes[temp[j]].energy -
						newEk + 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
				}
				else
				{
					double newEk = 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
					old_ext.mass += Area * fluxes[temp[j]].mass;
					old_ext.momentum += Area * fluxes[temp[j]].momentum;
					old_ext.energy += Area * fluxes[temp[j]].energy;
					old_ext.internal_energy += Area * fluxes[temp[j]].energy +
						newEk - 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
				}
			}
			for (size_t j = 0; j < temp.size(); ++j)
			{
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				double Area = tess.GetArea(temp[j]) * dt;
				double Ek = 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
				delta = Area * fluxes[temp[j]];
				double dEtherm = 0;
				if (N1 == i)
				{
					old_ext += delta;
					double Eknew = 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
					dEtherm = delta.energy - (Eknew - Ek);
					old_ext.internal_energy += delta.energy - (Eknew - Ek);
				}
				else
				{
					old_ext -= delta;
					double Eknew = 0.5 * ScalarProd(old_ext.momentum, old_ext.momentum) / old_ext.mass;
					dEtherm = -delta.energy - (Eknew - Ek);
					old_ext.internal_energy += dEtherm;
				}
				std::cout << "Face " << temp[j] << " neigh " << N0 << "," << N1 << " mass=" << fluxes[temp[j]].mass * Area <<
					" energy= " << fluxes[temp[j]].energy * Area << " Etherm= " << dEtherm << " momentum= " << abs(fluxes[temp[j]].momentum) * Area <<
					" Area*dt " << Area << " normal " << tess.Normal(temp[j]).x << "," << tess.Normal(temp[j]).y << "," << tess.Normal(temp[j]).z << std::endl;
				for (size_t k = 0; k < tracerstickernames.tracer_names.size(); ++k)
				{
					std::cout << tracerstickernames.tracer_names[k] << " flux is " << fluxes[temp[j]].tracers[k] * Area << std::endl;
				}
			}
			for (size_t j = 0; j < temp.size(); ++j)
			{
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				size_t Nother = N1 == i ? N0 : N1;
				std::cout << "Neigh cell " << Nother << ", density " << cells[Nother].density << " pressure " <<
					cells[Nother].pressure << " vx " << cells[Nother].velocity.x << " vy " << cells[Nother].velocity.y << " vz "
					<< cells[Nother].velocity.z << std::endl;
				for (size_t k = 0; k < tracerstickernames.tracer_names.size(); ++k)
				{
					std::cout << tracerstickernames.tracer_names[k] << " of other " << cells[Nother].tracers[k] << std::endl;
				}
			}
			assert(false);
		}
	}
	extensives.resize(tess.GetPointNo());
}


void RegularExtensiveUpdate3D::operator()(const vector<Conserved3D>& /*fluxes*/, const Tessellation3D& /*tess*/, const double /*dt*/,
	const vector<ComputationalCell3D>& /*cells*/, vector<Conserved3D>& /*extensives*/,
	size_t /*index*/, double /*time*/,
	TracerStickerNames const& /*tracerstickernames*/)const
{
	//assert(extensives[index].internal_energy > 0);
	return;
}
