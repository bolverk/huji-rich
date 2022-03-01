#include "Diffusion.hpp"

namespace
{
    // Vector3D calc_grad(double const cell,
	// 	Vector3D const& center, Vector3D const& cell_cm, double cell_volume, vector<double> const& neighbors,
	// 	vector<Vector3D> const& neighbor_centers, vector<Vector3D> const& neigh_cm, Tessellation3D const& tess,
	// 	Slope3D &temp, face_vec const& faces, std::vector<Vector3D>& c_ij)
	// {
    //     Vector3D v_temp(0, 0, 0);
	// 	size_t n = neighbor_centers.size();
	// 	// Create the matrix to invert and the vector to compare
	// 	std::array<double, 9>  m;
	// 	std::fill_n(m.begin(), 9, 0.0);
	// 	Vector3D c_ij;

	// 	for (size_t i = 0; i < n; i++)
	// 	{
	// 		c_ij = neigh_cm[i];
	// 		c_ij += cell_cm;
	// 		c_ij *= -0.5;
	// 		c_ij += tess.FaceCM(faces[i]);
	// 		const Vector3D r_ij = normalize(neighbor_centers[i] - center);
	// 		const double A = tess.GetArea(faces[i]);
	// 		m[0] -= c_ij.x*r_ij.x*A;
	// 		m[1] -= c_ij.y*r_ij.x*A;
	// 		m[2] -= c_ij.z*r_ij.x*A;
	// 		m[3] -= c_ij.x*r_ij.y*A;
	// 		m[4] -= c_ij.y*r_ij.y*A;
	// 		m[5] -= c_ij.z*r_ij.y*A;
	// 		m[6] -= c_ij.x*r_ij.z*A;
	// 		m[7] -= c_ij.y*r_ij.z*A;
	// 		m[8] -= c_ij.z*r_ij.z*A;
    //         v_temp += (A * (cell + neighbors[i])) * r_ij;
	// 	}
	// 	double const v_inv = 1.0 / cell_volume;
	// 	for (size_t i = 0; i < 9; ++i)
	// 		m[i] *= v_inv;
	// 	m[0] += 1;
	// 	m[4] += 1;
	// 	m[8] += 1;
	// 	// Find the det
	// 	const double det = -m[2] * m[4] * m[6] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[0] * m[5] * m[7] -
	// 		m[1] * m[3] * m[8] + m[0] * m[4] * m[8];
	// 	// Check none singular
	// 	if (std::abs(det) < 1e-10)
	// 	{
	// 		UniversalError eo("Singular matrix in calc_grad");
	// 		throw eo;
	// 	}
	// 	// Invert the matrix
	// 	std::array<double, 9>  m_inv;
	// 	std::fill_n(m_inv.begin(), 9, 0);
	// 	m_inv[0] = m[4] * m[8] - m[5] * m[7];
	// 	m_inv[1] = m[2] * m[7] - m[1] * m[8];
	// 	m_inv[2] = m[1] * m[5] - m[2] * m[4];
	// 	m_inv[3] = m[5] * m[6] - m[3] * m[8];
	// 	m_inv[4] = m[0] * m[8] - m[2] * m[6];
	// 	m_inv[5] = m[2] * m[3] - m[5] * m[0];
	// 	m_inv[6] = m[3] * m[7] - m[6] * m[4];
	// 	m_inv[7] = m[6] * m[1] - m[0] * m[7];
	// 	m_inv[8] = m[4] * m[0] - m[1] * m[3];
    //     double const det_volume_inverse = 1.0 / (2 * cell_volume * det);
	// 	for (size_t i = 0; i < 9; ++i)
	// 		m_inv[i] *= det_volume_inverse;
    //     Vector3D res(v_temp.x * m_inv[0] + v_temp.y * m_inv[1] + v_temp.z * m_inv[2],
    //         v_temp.x * m_inv[3] + v_temp.y * m_inv[4] + v_temp.z * m_inv[5],
    //         v_temp.x * m_inv[6] + v_temp.y * m_inv[7] + v_temp.z * m_inv[8]);
    //     return res;
	// }

    // Larsen second order flux limiter, taken from eq. 10 in "Diffusion, P1, and other approximate forms of radiation transport"
    double CalcSingleFluxLimiter(Vector3D const& grad, double const D, double const cell_value)
    {
        return 1.0 / (1 + ScalarProd(grad, grad) * 9 * D * D / (cell_value * cell_value * CG::speed_of_light * CG::speed_of_light));
    }
}

void Diffusion::BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells, std::string const& key_name,
    double const dt, std::vector<double>& b, std::vector<double>& x0) const
{
    size_t const key_index = binary_index_find(ComputationalCell3D::tracerNames, key_name);
    size_t const Nlocal = tess.GetPointNo();
    b.resize(Nlocal, 0);
    x0.resize(Nlocal, 0);
    std::vector<double> D(Nlocal), fleck_factor(Nlocal);
    std::vector<size_t> neighbors;
    face_vec faces;
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i);
        double const Er = cells[i].tracers[key_index] * cells[i].density;
        b[i] = Er * volume;
        x0[i] = Er;
        D[i] = D_coefficient_calcualtor_.CalcDiffusionCoefficient(i, cells);
        double const T = cells[i].temperature;
        double const sigma_planck = D_coefficient_calcualtor_.CalcPlanckOpacity(i, cells);
        double const beta = 4 * CG::radiation_constant * T * T * T / (cells[i].density * eos_.dT2cv(cells[i].density, T, cells[i].tracers, ComputationalCell3D::tracerNames));
        fleck_factor[i] = 1.0 / (1 +sigma_planck * CG::speed_of_light * dt * beta);
        b[i] += volume * fleck_factor[i] * dt * CG::speed_of_light * sigma_planck * T * T * T * T * CG::radiation_constant;
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
        double const sigma_planck = D_coefficient_calcualtor_.CalcPlanckOpacity(i, cells);
        double const T = cells[i].temperature;
        A[i].push_back(volume * (1 + fleck_factor[i] * dt * CG::speed_of_light * sigma_planck));
    }

    for(size_t i = 0; i < Nlocal; ++i)
    {
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const CM = tess.GetCellCM(i);
        Vector3D const point = tess.GetMeshPoint(i); 
        double const Dcell = D[i];
        double const Er = cells[i].tracers[key_index] * cells[i].density;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            // Here we assume no flux to outside cells, this needs to be changed to a general boundary condition
            size_t const neighbor_j = neighbors[j];
            if(i < neighbor_j)
            {
                if(!tess.IsPointOutsideBox(neighbor_j))
                {
                    double const Er_j = cells[neighbor_j].tracers[key_index] * cells[neighbor_j].density;
                    Vector3D const cm_ij = CM - tess.GetCellCM(neighbor_j);
                    Vector3D const grad_E = cm_ij * (1.0 / ScalarProd(cm_ij, cm_ij));
                    Vector3D const r_ij = point - tess.GetMeshPoint(neighbor_j);
                    double const flux_limiter = CalcSingleFluxLimiter(grad_E * (Er - Er_j), Dcell, 0.5 * (Er + Er_j));
                    double const mid_D = 0.5 * (D[neighbor_j] + Dcell) * flux_limiter;
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