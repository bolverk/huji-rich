#include "DiffusionForce.hpp"
#include <boost/math/special_functions/pow.hpp>
// equations taken from "EQUATIONS AND ALGORITHMS FOR MIXED-FRAME FLUX-LIMITED DIFFUSION RADIATION HYDRODYNAMICS"

void DiffusionForce::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes,const vector<Vector3D>& point_velocities, const double t,double dt,
		vector<Conserved3D> &extensives) const
{
    int total_iters = 0;
	double const CG_eps = 1e-8;
	std::vector<double> new_Er = CG::conj_grad_solver(CG_eps, total_iters, tess, cells , key_, dt, diffusion_);
	size_t const N = tess.GetPointNo();
	double max_Er = *std::max_element(new_Er.begin(), new_Er.end());
#ifdef RICH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &max_Er, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
	double max_diff = 0;
    size_t max_loc = 0;
    size_t const key_index = binary_index_find(ComputationalCell3D::tracerNames, key_);
	for(size_t i = 0; i < N; ++i)
	{
		double const diff = std::abs(new_Er[i] - cells[i].tracers[key_index] * cells[i].density) / (new_Er[i] + 0.01 * max_Er);
        if(diff > max_diff)
        {
            max_diff = diff;
            max_loc = i;
        }
	}
#ifdef RICH_MPI
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::pair<double, int> max_data(max_diff, rank);         
    MPI_Allreduce(MPI_IN_PLACE, &max_data, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#endif
    std::vector<size_t> neighbors;
    face_vec faces;
    std::vector<double> flux_limiter(N, 0), R2(N, 0);
    for(size_t i = 0; i < N; ++i)
    {
        double const volume = tess.GetVolume(i);
        // The transport contribution
        extensives[i].tracers[key_index] = new_Er[i] * volume;
        double const T = cells[i].temperature;
        // The contribution from matter-radiation coupling
        double const dE = diffusion_.fleck_factor[i] * CG::speed_of_light * dt * diffusion_.sigma_planck[i] * (new_Er[i] - T * T * T * T * CG::radiation_constant) * volume;
        extensives[i].energy += dE;
        extensives[i].internal_energy += dE;
        // Calcualte gradient of radiation field
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const point = tess.GetMeshPoint(i);
        Vector3D gradE(0, 0, 0);
        for(size_t j = 0; j < Nneigh; ++j)
        {
            size_t const neighbor_j = neighbors[j];
            Vector3D r_ij = point - tess.GetMeshPoint(neighbor_j);
            r_ij *= 1.0 / abs(r_ij);
            double Emid = 0;           
            if(!tess.IsPointOutsideBox(neighbor_j))
                Emid = 0.5 * (new_Er[i] + new_Er[neighbor_j]);
            else
            {
                Vector3D dummy_v;
                diffusion_.boundary_calc_.GetOutSideValues(tess, cells, i, neighbor_j, new_Er, Emid, dummy_v, key_index);
                Emid *= 0.5;
                Emid += 0.5 * new_Er[i];
            }
            gradE += r_ij * (tess.GetArea(faces[j]) * Emid);
        }
        gradE *= 1.0 / volume;
        double const D = diffusion_.D_coefficient_calcualtor.CalcDiffusionCoefficient(cells[i]);
        flux_limiter[i] = CG::CalcSingleFluxLimiter(gradE, D, new_Er[i]);
        R2[i] = flux_limiter[i] / 3+ boost::math::pow<2>(flux_limiter[i] * abs(gradE) * D 
            / (CG::speed_of_light * new_Er[i]));
        // Add radition pressure term to momentum
        double const old_Ek = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
        extensives[i].momentum -= (volume * flux_limiter[i] * dt / 3) * gradE;
        double const new_Ek = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
        // Add work done by radiation on gas and vice versa
        extensives[i].energy += new_Ek - old_Ek;
        extensives[i].tracers[key_index] -= new_Ek - old_Ek;
        // Add relativity term
        double const dE_rad = dt * flux_limiter[i] * 2 * diffusion_.sigma_planck[i] * D / CG::speed_of_light * ScalarProd(cells[i].velocity, gradE) * volume;
        extensives[i].internal_energy += dE_rad;
        extensives[i].energy += dE_rad;
        extensives[i].tracers[key_index] -= dE_rad;
    }
#ifdef RICH_MPI
    MPI_exchange_data2(tess, R2, true);
#endif
    for(size_t i = 0; i < N; ++i)
    {
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const point = tess.GetMeshPoint(i);
        double dE = 0;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            size_t const neighbor_j = neighbors[j];
            Vector3D r_ij = point - tess.GetMeshPoint(neighbor_j);
            r_ij *= 1.0 / abs(r_ij);
            Vector3D velocity_outside;
            double Er_outside, R2_outside;
            // Add enthalpy advection, remember that we already had some advection in the hydro
            if(!tess.IsPointOutsideBox(neighbor_j))
            {
                velocity_outside = cells[neighbor_j].velocity;
                Er_outside = new_Er[neighbor_j];
                R2_outside = R2[neighbor_j];
            }
            else
            {
                diffusion_.boundary_calc_.GetOutSideValues(tess, cells, i, neighbor_j, new_Er, Er_outside, velocity_outside, key_index);
                R2_outside = R2[i];
            }
            Vector3D const velocity_mid = 0.5 * (cells[i].velocity + velocity_outside);
            double const v_normal = ScalarProd(velocity_mid, r_ij);
            if(v_normal > 0)
                dE += Er_outside * tess.GetArea(faces[j]) * dt * v_normal * (0.5 - 0.5 * R2_outside);
            else
                dE += new_Er[i] * tess.GetArea(faces[j]) * dt * v_normal * (0.5 - 0.5 * R2[i]);
        }
        extensives[i].tracers[key_index] += dE ;

        dE = (1.5 - 0.5 * R2[i]) * diffusion_.sigma_planck[i] * tess.GetVolume(i) * dt *
            new_Er[i] * ScalarProd(cells[i].velocity, cells[i].velocity) / CG::speed_of_light;
        extensives[i].tracers[key_index] += dE ;
        extensives[i].energy += dE;
        extensives[i].internal_energy += dE;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy) || extensives[i].tracers[key_index] < 0)
            throw UniversalError("Bad internal energy in DiffusionForce::operator()");
    }
    diffusion_.fleck_factor.clear();
    diffusion_.fleck_factor.shrink_to_fit();
    diffusion_.sigma_planck.clear();
    diffusion_.sigma_planck.shrink_to_fit();
	next_dt_ = dt * std::min(0.1 / max_diff, 1.1);
}   


double DiffusionForce::SuggestInverseTimeStep(void)const
{
    return 1.0 / next_dt_;
}