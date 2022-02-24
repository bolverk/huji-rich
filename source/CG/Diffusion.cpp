#include "Diffusion.hpp"

namespace
{
    double CalcFluxLimiter(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const key_index, double const D,
        std::vector<size_t> const& neighbors, face_vec const& faces)
    {
        // Calculate grad(key)
        Vector3D const point = tess.GetMeshPoint(index); 
        size_t const Nneigh = neighbors.size();
        double const cell_value = cells[index].tracers[key_index] * cells[index].density;
        Vector3D grad_key;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            // Here we ignore outside cells, this needs to be changed to a general boundary condition
            size_t const neighbor_j = neighbors[j];
            Vector3D const r_ij = point - tess.GetMeshPoint(neighbor_j);
            if(!tess.IsPointOutsideBox(neighbor_j))
                grad_key += r_ij * (tess.GetArea(faces[j]) * 0.5 * (cell_value + cells[neighbor_j].tracers[key_index] * cells[neighbor_j].density) / abs(r_ij));
            else
                grad_key += r_ij * (tess.GetArea(faces[j]) * cell_value  / abs(r_ij));
        }
        grad_key *= 1.0 / tess.GetVolume(index);
        return 1.0 / std::sqrt(1 + D * D * ScalarProd(grad_key, grad_key) / (cell_value * cell_value * CG::speed_of_light * CG::speed_of_light));
    }
}

void Diffusion::BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells, std::string const& key_name,
    double const dt, std::vector<double>& b, std::vector<double>& x0) const
{
    size_t const key_index = binary_index_find(ComputationalCell3D::tracerNames, key_name);
    size_t const Nlocal = tess.GetPointNo();
    b.resize(Nlocal, 0);
    x0.resize(Nlocal, 0);
    std::vector<double> D(Nlocal), flux_limiter(Nlocal);
    std::vector<size_t> neighbors;
    face_vec faces;
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i);
        double const Er = cells[i].tracers[key_index] * cells[i].density;
        b[i] = Er * volume;
        x0[i] = Er;
        D[i] = D_coefficient_calcualtor_.CalcDiffusionCoefficient(i, cells);
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        flux_limiter[i] = CalcFluxLimiter(tess, cells, i, key_index, D[i], neighbors, faces);
        D[i] *= flux_limiter[i];
    }
#ifdef RICH_MPI
    MPI_exchange_data2(tess, D, true);
#endif
    size_t max_neigh = 0;
    // Find maximum number of neighbors and allocate data
    for(size_t i = 0; i < Nlocal; ++i)
        max_neigh = std::max(max_neigh, tess.GetNeighbors(i).size());
    ++max_neigh;
    A.clear();
    A.resize(Nlocal);
    A_indeces.clear();
    A_indeces.resize(Nlocal);

    // Build the matrix
    for(size_t i = 0; i < Nlocal; ++i)
    {
        A_indeces[i].push_back(i);
        double const volume = tess.GetVolume(i);
        A[i].push_back(volume);
    }

    for(size_t i = 0; i < Nlocal; ++i)
    {
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const CM = tess.GetCellCM(i);
        Vector3D const point = tess.GetMeshPoint(i); 
        double const Dcell = D[i];
        for(size_t j = 0; j < Nneigh; ++j)
        {
            // Here we assume no flux to outside cells, this needs to be changed to a general boundary condition
            size_t const neighbor_j = neighbors[j];
            if(i < neighbor_j)
            {
                if(!tess.IsPointOutsideBox(neighbor_j))
                {
                    Vector3D const cm_ij = CM - tess.GetCellCM(neighbor_j);
                    Vector3D const grad_E = cm_ij * (1.0 / ScalarProd(cm_ij, cm_ij));
                    Vector3D const r_ij = point - tess.GetMeshPoint(neighbor_j);
                    double const mid_D = 0.5 * (D[neighbor_j] + Dcell);
                    double const flux = ScalarProd(grad_E, r_ij * (tess.GetArea(faces[j]) / abs(r_ij))) * dt * mid_D; 
                    if(neighbor_j < Nlocal)
                    {
                        A[i][0] += flux;
                        A[i].push_back(-flux);
                        A_indeces[i].push_back(neighbor_j);
                        A[neighbor_j].push_back(-flux);
                        A_indeces[neighbor_j].push_back(i);
                        A[neighbor_j][0] += flux;
                    }                  
                    else
                    {
                        A[i][0] += flux;
                        A[i].push_back(-flux);
                        A_indeces[i].push_back(neighbor_j);
                    }
                }
                else
                    boundary_calc_.SetBoundaryValues(tess, i, neighbor_j, dt, cells, key_index, tess.GetArea(faces[j]), A[i][0], b[i]);
            }
        }
    }
    for(size_t i = 0; i < Nlocal; ++i)
    {
        A[i].resize(max_neigh, 0);
        A_indeces[i].resize(max_neigh, max_size_t);
    }
}

void DiffusionSideBoundary::SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt, 
        std::vector<ComputationalCell3D> const& /*cells*/, size_t const /*key_index*/, double const Area, double& A, double &b)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        A += 0.5 * CG::speed_of_light * dt * Area;
        b += 2 * Area * dt * CG::stefan_boltzman * T_ * T_ * T_ * T_;
    }
}

double PowerLawOpacity::CalcDiffusionCoefficient(size_t const index, std::vector<ComputationalCell3D> const& cells) const
{
    double const T = std::pow(cells[index].tracers[0] * cells[index].density / CG::radiation_constant, 0.25);
    return D0_ * std::pow(cells[index].density, alpha_) * std::pow(T, beta_);
}

double PowerLawOpacity::CalcPlanckOpacity(size_t const index, std::vector<ComputationalCell3D> const& cells) const
{
    return CalcDiffusionCoefficient(index, cells);
}