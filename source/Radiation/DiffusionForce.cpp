#include "DiffusionForce.hpp"
#include <boost/math/special_functions/pow.hpp>
// equations taken from "EQUATIONS AND ALGORITHMS FOR MIXED-FRAME FLUX-LIMITED DIFFUSION RADIATION HYDRODYNAMICS"

void DiffusionForce::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes,const vector<Vector3D>& point_velocities, const double t,double dt,
		vector<Conserved3D> &extensives) const
{
 #ifdef RICH_MPI
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    std::vector<Conserved3D> old_extensives(extensives);
    std::vector<size_t> neighbors;
    face_vec faces;
	size_t const N = tess.GetPointNo();
    std::vector<double> flux_limiter(N, 0), R2(N, 0);
    std::vector<double> new_Er(N, 0);
    for(size_t i = 0; i < N; ++i)
        new_Er[i] = cells[i].Erad * cells[i].density;
    if(N == 0)
        std::cout<<"Zero nubmer of cells in DiffForce"<<std::endl;
	double max_Er = *std::max_element(new_Er.begin(), new_Er.end());
#ifdef RICH_MPI
    MPI_exchange_data2(tess, new_Er, true);
    MPI_Allreduce(MPI_IN_PLACE, &max_Er, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    size_t const Nzero = diffusion_.zero_cells_.size();
    std::vector<size_t> zero_indeces;
    for(size_t i = 0; i < Nzero; ++i)
        zero_indeces.push_back(binary_index_find(ComputationalCell3D::stickerNames, diffusion_.zero_cells_[i]));
    double min_dt_inv = std::numeric_limits<double>::min() * 100;
    size_t min_dt_inv_index = 0;

    for(size_t i = 0; i < N; ++i)
    {
        bool to_calc = true;
        for(size_t j = 0; j < Nzero; ++j)
            if(cells[i].stickers[zero_indeces[j]])
                to_calc = false;
        if(not to_calc)
            continue;
        double const volume = tess.GetVolume(i);
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
                diffusion_.boundary_calc_.GetOutSideValues(tess, cells, i, neighbor_j, new_Er, Emid, dummy_v);
                Emid *= 0.5;
                Emid += 0.5 * new_Er[i];
            }
            gradE += r_ij * (tess.GetArea(faces[j]) * Emid);
        }
        gradE *= -1.0 / volume;
        double const D = diffusion_.D_coefficient_calcualtor.CalcDiffusionCoefficient(cells[i]);
        flux_limiter[i] = diffusion_.flux_limiter_ ? CG::CalcSingleFluxLimiter(gradE, D, new_Er[i]) : 1;
        R2[i] = diffusion_.flux_limiter_ ? flux_limiter[i] / 3 + boost::math::pow<2>(flux_limiter[i] * abs(gradE) * D 
            / (CG::speed_of_light * new_Er[i])) : 1.0 / 3.0;
        if(not momentum_limit_)
        {
            flux_limiter[i] = 1;
            R2[i] = 1.0 / 3.0;
        }
        // Add radition pressure term to momentum
        double const old_Ek = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
        extensives[i].momentum -= (volume * flux_limiter[i] * dt / 3) * gradE;
        double const new_Ek = 0.5 * ScalarProd(extensives[i].momentum, extensives[i].momentum) / extensives[i].mass;
        // Add work done by radiation on gas and vice versa
        extensives[i].energy += new_Ek - old_Ek;
        extensives[i].Erad -= new_Ek - old_Ek;
        // Add relativity term
        double const dE_rad = std::min(1.0 / 3.0, diffusion_.D_coefficient_calcualtor.CalcPlanckOpacity(cells[i]) * D / CG::speed_of_light) * dt * flux_limiter[i] * 2 * ScalarProd(cells[i].velocity, gradE) * volume;
        double const new_min_dt_inv = std::abs(dE_rad / extensives[i].Erad);
        if(new_min_dt_inv > min_dt_inv)
        {
            min_dt_inv = new_min_dt_inv;
            min_dt_inv_index = i;
        }
        extensives[i].internal_energy += dE_rad;
        extensives[i].energy += dE_rad;
        extensives[i].Erad -= dE_rad;
        //protect against negative Erad
        if(not diffusion_.flux_limiter_ && extensives[i].Erad < 0 && new_Er[i] < max_Er * 1e-6)
            extensives[i].Erad = old_extensives[i].Erad;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy) || extensives[i].Erad < 0)
        {
            UniversalError eo("Negative energy in DiffusionForce");
            eo.addEntry("new_Er", new_Er[i]);
            eo.addEntry("internal_energy", extensives[i].internal_energy);
            eo.addEntry("Erad", extensives[i].Erad);
            eo.addEntry("X", point.x);
            throw eo;
        }
    }
#ifdef RICH_MPI
    MPI_exchange_data2(tess, R2, true);
#endif
    for(size_t i = 0; i < N; ++i)
    {
        bool to_calc = true;
        for(size_t j = 0; j < Nzero; ++j)
            if(cells[i].stickers[zero_indeces[j]])
                to_calc = false;
        if(not to_calc)
            continue;
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
            double Er_outside, R2_outside, density_outside;
            // Add enthalpy advection, remember that we already had some advection in the hydro
            if(!tess.IsPointOutsideBox(neighbor_j))
            {
                velocity_outside = cells[neighbor_j].velocity;
                Er_outside = new_Er[neighbor_j];
                R2_outside = R2[neighbor_j];
                density_outside = cells[neighbor_j].density;
            }
            else
            {
                diffusion_.boundary_calc_.GetOutSideValues(tess, cells, i, neighbor_j, new_Er, Er_outside, velocity_outside);
                R2_outside = R2[i];
                density_outside = cells[i].density;
            }
            size_t n0 = tess.GetFaceNeighbors(faces[j]).first;
            double dm = fluxes[faces[j]].mass;
            if(i == n0)
                dm *= -1;
            if(dm > 0)
                dE += Er_outside * tess.GetArea(faces[j]) * dt * dm * (0.5 - 0.5 * R2_outside) / density_outside;
            else
                dE += new_Er[i] * tess.GetArea(faces[j]) * dt * dm * (0.5 - 0.5 * R2[i]) / cells[i].density;
        }
        extensives[i].Erad += dE ;
        if(not diffusion_.flux_limiter_ && extensives[i].Erad < 0 && new_Er[i] < max_Er * 1e-6)
            extensives[i].Erad = old_extensives[i].Erad;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy) || extensives[i].Erad < 0)
            throw UniversalError("Negative energy in DiffusionForce2");
        double const D = diffusion_.D_coefficient_calcualtor.CalcDiffusionCoefficient(cells[i]);
        dE = (1.5 - 0.5 * R2[i]) * std::min(diffusion_.D_coefficient_calcualtor.CalcPlanckOpacity(cells[i]), 3 * CG::speed_of_light / D) * tess.GetVolume(i) * dt *
            new_Er[i] * ScalarProd(cells[i].velocity, cells[i].velocity) / CG::speed_of_light;
        extensives[i].Erad += dE ;
        extensives[i].energy -= dE;
        extensives[i].internal_energy -= dE;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy) || extensives[i].Erad < 0)
            throw UniversalError("Bad internal energy in DiffusionForce::operator()");
    }


    double max_diff = 0;
    size_t max_loc = 0;
    for(size_t i = 0; i < N; ++i)
    {
        double const diff = std::abs(new_Er[i] - cells[i].Erad * cells[i].density) / (new_Er[i] + 0.005 * max_Er);
        if(diff > max_diff)
        {
            max_diff = diff;
            max_loc = i;
        }
    }
#ifdef RICH_MPI
    std::pair<double, int> max_data(max_diff, rank);         
    MPI_Allreduce(MPI_IN_PLACE, &max_data, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    max_diff = max_data.first;
#endif
	next_dt_ = dt * std::min(0.05 / max_diff, 1.1);
}   


double DiffusionForce::SuggestInverseTimeStep(void)const
{
    return 1.0 / next_dt_;
}