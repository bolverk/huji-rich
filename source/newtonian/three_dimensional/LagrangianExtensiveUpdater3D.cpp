#include "LagrangianExtensiveUpdater3D.hpp"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

LagrangianExtensiveUpdater3D::LagrangianExtensiveUpdater3D(LagrangianFlux3D const & lflux, EquationOfState const & eos,
	Ghost3D const & ghost, const vector<pair<const ConditionExtensiveUpdater3D::Condition3D*,
	const ConditionExtensiveUpdater3D::Action3D*>>&sequence) :lflux_(lflux), eos_(eos), ghost_(ghost), sequence_(sequence) {}

void LagrangianExtensiveUpdater3D::operator()(const vector<Conserved3D>& fluxes, const Tessellation3D & tess,
	const double dt, const vector<ComputationalCell3D>& cells, vector<Conserved3D>& extensives, double time,
					      TracerStickerNames const & tracerstickernames, std::vector<Vector3D> const& /*face_vel*/, 
					      std::vector<std::pair<ComputationalCell3D, ComputationalCell3D> > const& /*interp_values*/) const
{
	std::vector<Conserved3D> old_extensive(extensives);
#ifdef RICH_MPI
	Conserved3D edummy;
	MPI_exchange_data(tess, old_extensive, true,&edummy);
#endif
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	size_t indexZ = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaZ")) - tracerstickernames.tracer_names.begin());

	size_t N = tess.GetPointNo();
	// Reduce the area tracer
	bool AreaChange = indexX < tracerstickernames.tracer_names.size();
	if (AreaChange)
	{
		for (size_t i = 0; i < N; ++i)
		{
			extensives[i].tracers[indexX] *= 0.95;
			extensives[i].tracers[indexY] *= 0.95;
			extensives[i].tracers[indexZ] *= 0.95;
		}
	}

	std::vector<double> oldEk(N, 0), oldEtherm(N, 0), oldE(N, 0);
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
		double FArea = dt*tess.GetArea(i);
		delta = fluxes[i] * FArea;
		size_t n0 = tess.GetFaceNeighbors(i).first;
		size_t n1 = tess.GetFaceNeighbors(i).second;

		double deltaWs = -lflux_.ws_[i] * tess.GetArea(i) * dt;
		Vector3D normal = normalize(tess.GetMeshPoint(n1) - tess.GetMeshPoint(n0));
		if (lflux_.Lag_calc_[i])
		{
			double p_star = ScalarProd(fluxes[i].momentum, normal);
			if (p_star > 0)
			{
				double v_star = fluxes[i].energy / p_star;
				double v_new = (v_star - lflux_.ws_[i]);
				if (v_new*v_star > 0)
				{
					if (v_new > 0 && !tess.IsPointOutsideBox(n1) && p_star > 1.2*cells[n1].pressure)
					{
						double Ek = 0.5*ScalarProd(old_extensive[n1].momentum, old_extensive[n1].momentum) / old_extensive[n1].mass;
						double Eknew = 0.5*ScalarProd(old_extensive[n1].momentum + delta.momentum,
							old_extensive[n1].momentum + delta.momentum) / old_extensive[n1].mass;
						v_new = std::max(v_new, (Eknew - Ek) / std::max(1e-50, p_star*FArea));
					}
					if (v_new < 0 && !tess.IsPointOutsideBox(n0) && p_star > 1.2*cells[n0].pressure)
					{
						double Ek = 0.5*ScalarProd(old_extensive[n0].momentum, old_extensive[n0].momentum) /
							old_extensive[n0].mass;
						double Eknew = 0.5*ScalarProd(old_extensive[n0].momentum - delta.momentum,
							old_extensive[n0].momentum - delta.momentum) / old_extensive[n0].mass;
						v_new = std::min(v_new, (Ek - Eknew) / std::max(1e-50, p_star*FArea));
					}
					delta.energy = p_star*v_new*tess.GetArea(i) * dt;
				}
				else
				{
					if (v_new > 0 && !tess.IsPointOutsideBox(n0))
					{
						double newp = std::min(p_star, cells[n0].pressure);
						delta.energy = tess.GetArea(i) * dt*v_new*newp;
						delta.momentum *= newp / p_star;
					}
					else
						if (!tess.IsPointOutsideBox(n1))
						{
							double newp = std::min(p_star, cells[n1].pressure);
							delta.energy = tess.GetArea(i) * dt*v_new*newp;
							delta.momentum *= newp / p_star;
						}
						else
							delta.energy = 0;
				}
			}
		}
		if (n0 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0] -= delta;
			double Eknew = 0.5*ScalarProd(extensives[n0].momentum, extensives[n0].momentum) / extensives[n0].mass;
			extensives[n0].internal_energy -= delta.energy + (Eknew - Ek);
			if (AreaChange)
			{
				extensives[n0].tracers[indexX] -= normal.x*deltaWs;
				extensives[n0].tracers[indexY] -= normal.y*deltaWs;
				extensives[n0].tracers[indexZ] -= normal.z*deltaWs;
			}
		}
		if (n1 < N)
		{
			double Ek = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1] += delta;
			double Eknew = 0.5*ScalarProd(extensives[n1].momentum, extensives[n1].momentum) / extensives[n1].mass;
			extensives[n1].internal_energy += delta.energy - (Eknew - Ek);
			if (AreaChange)
			{
				extensives[n1].tracers[indexX] -= normal.x*deltaWs;
				extensives[n1].tracers[indexY] -= normal.y*deltaWs;
				extensives[n1].tracers[indexZ] -= normal.z*deltaWs;
			}
		}
	}

	for (size_t i = 0; i < N; ++i)
	{
		double dEtherm = extensives[i].internal_energy - oldEtherm[i];
		double dEk = 0.5*ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass - oldEk[i];
		double dE = extensives[i].energy - oldE[i];
		if (dEtherm*(dE - dEk) > 0)
			if (std::abs(dEtherm) > 0.999 *std::abs(dE - dEk) && std::abs(dEtherm) < 1.001*std::abs(dE - dEk))
				extensives[i].internal_energy = oldEtherm[i] + (dE - dEk);
		// check cell
		if ((!(extensives[i].mass > 0)) || (!(extensives[i].energy > 0)) || (!std::isfinite(extensives[i].momentum.x)) || (!std::isfinite(extensives[i].momentum.y))
			|| (!std::isfinite(extensives[i].momentum.z)))
		{
			int rank = 0;
#ifdef RICH_MPI
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
			std::cout << "Bad cell in LagrangianExtensiveUpdate, cell " << i << " rank " << rank <<" time "<<time<< std::endl;
			std::cout << "mass " << extensives[i].mass << " energy " << extensives[i].energy << " internalE " <<
				extensives[i].internal_energy << " momentum" << abs(extensives[i].momentum) << " volume " << tess.GetVolume(i)
				<< std::endl;
			std::cout << "old mass " << old_extensive[i].mass << " energy " << old_extensive[i].energy << " internalE " <<
				old_extensive[i].internal_energy << " momentum" << abs(old_extensive[i].momentum) << " volume " << tess.GetVolume(i)
				<< std::endl;
			std::cout << "Old cell, density " << cells[i].density << " pressure " << cells[i].pressure << " v " <<
				abs(cells[i].velocity) << std::endl;
			face_vec temp = tess.GetCellFaces(i);
			for (size_t j = 0; j < temp.size(); ++j)
			{
				size_t N0 = tess.GetFaceNeighbors(temp[j]).first;
				size_t N1 = tess.GetFaceNeighbors(temp[j]).second;
				double Area = tess.GetArea(temp[j]) * dt;
				std::cout << "Face " << temp[j] << " neigh " << N0 << "," << N1 << " mass=" << fluxes[temp[j]].mass*Area <<
					" energy " << fluxes[temp[j]].energy*Area << " momentum=" << abs(fluxes[temp[j]].momentum)*Area <<
					" Area*dt " << Area << " N0 "<<tess.GetMeshPoint(N0).x<<","<<tess.GetMeshPoint(N0).y<<","
					<<tess.GetMeshPoint(N0).z<<" N1 "<<tess.GetMeshPoint(N1).x<<","<<tess.GetMeshPoint(N1).y<<","
					<<tess.GetMeshPoint(N1).z<<" IDs "<<cells[N0].ID<<" "<<cells[N1].ID<< std::endl;
				std::cout << "dl " << cells[N0].density << " pl " << cells[N0].pressure << " vxl " << cells[N0].velocity.x << " vyl " << cells[N0].velocity.y << " vzl " << cells[N0].velocity.z << std::endl;
				std::cout << "dr " << cells[N1].density << " pr " << cells[N1].pressure << " vxr " << cells[N1].velocity.x << " vyr " << cells[N1].velocity.y << " vzr " << cells[N1].velocity.z << std::endl;
				Vector3D normal = normalize(tess.GetMeshPoint(N1) - tess.GetMeshPoint(N0));
				std::cout << "dx " << normal.x << " dy " << normal.y << " dz " << normal.z << std::endl;
				if (lflux_.Lag_calc_[temp[j]])
				{
					
					double p_star = ScalarProd(fluxes[temp[j]].momentum,normal);
					double v_star = fluxes[temp[j]].energy / p_star;
					double v_new = (v_star - lflux_.ws_[temp[j]]);
					std::cout << "Old pstar " << p_star << " vstar " << v_star <<" vnew "<<v_new<<std::endl;
					if (v_new*v_star > 0)
					{
						if (v_new > 0 && !tess.IsPointOutsideBox(N1) && p_star > 1.2*cells[N1].pressure)
						{
							double Ek = 0.5*ScalarProd(old_extensive[N1].momentum, old_extensive[N1].momentum) / 
								old_extensive[N1].mass;
							double Eknew = 0.5*ScalarProd(old_extensive[N1].momentum + delta.momentum,
								old_extensive[N1].momentum + delta.momentum) / old_extensive[N1].mass;
							v_new = std::max(v_new, (Eknew - Ek) / std::max(1e-50, p_star*Area));
						}
						if (v_new < 0 && !tess.IsPointOutsideBox(N0) && p_star > 1.2*cells[N0].pressure)
						{
							double Ek = 0.5*ScalarProd(old_extensive[N0].momentum, old_extensive[N0].momentum) / 
								old_extensive[N0].mass;
							double Eknew = 0.5*ScalarProd(old_extensive[N0].momentum - delta.momentum,
								old_extensive[N0].momentum - delta.momentum) / old_extensive[N0].mass;
							v_new = std::min(v_new, (Ek - Eknew) / std::max(1e-50, p_star*Area));
						}
					}
					else
					{
						if (v_new > 0 && !tess.IsPointOutsideBox(N0))
						{
							double newp = std::min(p_star, cells[N0].pressure);
							p_star = newp;
						}
						else
							if (!tess.IsPointOutsideBox(N1))
							{
								double newp = std::min(p_star, cells[N1].pressure);
								p_star = newp;
							}
							else
								v_new = 0;
					}
					std::cout << "Pstar " << p_star << " Ustar " << v_new << std::endl;
				}
			}
			assert(false);
		}
		for (size_t j = 0; j < sequence_.size(); ++j)
		{
			if (sequence_[j].first->operator()(i, tess, cells, time, tracerstickernames))
			{
				sequence_[j].second->operator()(fluxes, tess, dt, cells, extensives, i, time, tracerstickernames);
				break;
			}
		}
	}
	extensives.resize(tess.GetPointNo());
}
