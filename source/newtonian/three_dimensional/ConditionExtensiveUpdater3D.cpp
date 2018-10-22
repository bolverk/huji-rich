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

void ConditionExtensiveUpdater3D::operator()(const vector<Conserved3D>& fluxes,	const Tessellation3D& tess,
	const double dt,const vector<ComputationalCell3D>& cells,vector<Conserved3D>& extensives,double time,
	TracerStickerNames const& tracerstickernames) const
{
	size_t N = tess.GetPointNo();
	std::vector<double> oldEk(N, 0), oldEtherm(N, 0),oldE(N,0);
	for (size_t i = 0; i < N; ++i)
	{
		oldEk[i] = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
		oldEtherm[i] = extensives[i].internal_energy;
		oldE[i] = extensives[i].energy;
	}
	size_t Nfluxes = fluxes.size();
	Conserved3D delta;
	for (size_t i = 0; i < Nfluxes; ++i)
	{
		delta = fluxes[i] * dt*tess.GetArea(i);
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;
		if (n0 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0] -= delta;
			double Eknew = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0].internal_energy -= delta.energy + (Eknew - Ek);
		}
		if (n1 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1] += delta;
			double Eknew = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1].internal_energy += delta.energy - (Eknew - Ek);
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		double dEtherm = extensives[i].internal_energy - oldEtherm[i];
		double dEk = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass - oldEk[i];
		double dE = extensives[i].energy - oldE[i];
		if(dEtherm*(dE-dEk)>0)
			if (std::abs(dEtherm) > 0.999 *std::abs(dE - dEk) && std::abs(dEtherm) < 1.001*std::abs(dE - dEk))
				extensives[i].internal_energy = oldEtherm[i] + (dE - dEk);
		for (size_t j = 0; j < sequence_.size(); ++j)
		{
			if (sequence_[j].first->operator()(i, tess, cells, time, tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, tess, dt, cells, extensives, i, time, tracerstickernames);
				break;
			}
		}
		// check cell
		if (extensives[i].mass < 0 || extensives[i].energy < 0 || extensives[i].internal_energy < 0)
		{
			int rank = 0;
#ifdef RICH_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
			std::cout << "Bad cell in ExtensiveUpdater3D, cell " << i << " rank " << rank << std::endl;
			std::cout << "mass " << extensives[i].mass << " energy " << extensives[i].energy << " internalE " <<
				extensives[i].internal_energy << " momentum" << abs(extensives[i].momentum) << " volume " << tess.GetVolume(i)
				<< std::endl;
			std::cout << "Old cell, density " << cells[i].density << " pressure " << cells[i].pressure << " v " <<
				abs(cells[i].velocity) << std::endl;
			vector<size_t> temp = tess.GetCellFaces(i);
			for (size_t j = 0; j < temp.size(); ++j)
			{
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				double Area = tess.GetArea(temp[j]) * dt;
				std::cout << "Face " << temp[j] << " neigh " << N0 << "," << N1 << " mass=" << fluxes[temp[j]].mass*Area <<
					" energy " << fluxes[temp[j]].energy*Area << " momentum=" << abs(fluxes[temp[j]].momentum)*Area <<
					" Area*dt " << Area << std::endl;
			}
			for (size_t j = 0; j < temp.size(); ++j)
			{
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				size_t Nother = N1 == i ? N0 : N1;
				std::cout << "Neigh cell " <<Nother<<", density " << cells[Nother].density << " pressure " << 
					cells[Nother].pressure << " v " <<	abs(cells[Nother].velocity) << std::endl;
			}
			assert(false);
		}
	}
	extensives.resize(tess.GetPointNo());
}


void RegularExtensiveUpdate3D::operator()(const vector<Conserved3D>& /*fluxes*/, const Tessellation3D& /*tess*/, const double /*dt*/,
	const vector<ComputationalCell3D>& /*cells*/, vector<Conserved3D> &extensives, size_t index, double /*time*/,
	TracerStickerNames const& /*tracerstickernames*/)const
{
	assert(extensives[index].internal_energy > 0);
	return;
}