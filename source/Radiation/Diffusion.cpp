#include "Diffusion.hpp"

namespace CG
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
        double const R = std::max(3 * abs(grad) * D / (cell_value * CG::speed_of_light), 1e20 * std::numeric_limits<double>::min());
        if(R < 1e-2) //series expansion
            return 1 - R * R / 15 + 2 * R * R * R * R /315;
        return 3 * (1.0 / std::tanh(R) - 1.0 / R) / R;
    }
}

void Diffusion::BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells,
    double const dt, std::vector<double>& b, std::vector<double>& x0, double const current_time) const
{
    size_t const Nlocal = tess.GetPointNo();
    b.resize(Nlocal, 0);
    x0.resize(Nlocal, 0);
    D.resize(Nlocal);
    fleck_factor.resize(Nlocal);
    sigma_planck.resize(Nlocal);
    std::vector<size_t> neighbors;
    face_vec faces;
    std::vector<size_t> zero_indeces;
    size_t const Nzero = zero_cells_.size();
    for(size_t i = 0; i < Nzero; ++i)
        zero_indeces.push_back(binary_index_find(ComputationalCell3D::stickerNames, zero_cells_[i]));
    double const zero_value = 1e-10;
    for(size_t i = 0; i < Nlocal; ++i)
    {
        double const volume = tess.GetVolume(i);
        D[i] = D_coefficient_calcualtor.CalcDiffusionCoefficient(cells[i]);
        double const T = cells[i].temperature;
        sigma_planck[i] = D_coefficient_calcualtor.CalcPlanckOpacity(cells[i]);
        double const beta = 4 * CG::radiation_constant * T * T * T / (cells[i].density * eos_.dT2cv(cells[i].density, T, cells[i].tracers, ComputationalCell3D::tracerNames));
        fleck_factor[i] = 1.0 / (1 + sigma_planck[i] * CG::speed_of_light * dt * beta);
        bool set_to_zero = false;
        for(size_t j = 0; j < Nzero; ++j)
            if(cells[i].stickers[zero_indeces[j]])
                set_to_zero = true;
        double const Er = cells[i].Erad * cells[i].density * (set_to_zero ? zero_value : 1);
        b[i] = Er * volume;
        x0[i] = Er;       
        b[i] += volume * fleck_factor[i] * dt * CG::speed_of_light * sigma_planck[i] * T * T * T * T * CG::radiation_constant;
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
        double const T = cells[i].temperature;
        A[i].push_back(volume * (1 + fleck_factor[i] * dt * CG::speed_of_light * sigma_planck[i]));
    }

    for(size_t i = 0; i < Nlocal; ++i)
    {
        faces = tess.GetCellFaces(i);
        tess.GetNeighbors(i, neighbors);
        size_t const Nneigh = neighbors.size();
        Vector3D const CM = tess.GetCellCM(i);
        Vector3D const point = tess.GetMeshPoint(i); 
        double const Dcell = D[i];
        double const Er = cells[i].Erad * cells[i].density;
        bool self_zero = false;
        for(size_t k = 0; k < Nzero; ++k)
            if(cells[i].stickers[zero_indeces[k]])
                self_zero = true;
        for(size_t j = 0; j < Nneigh; ++j)
        {
            // Here we assume no flux to outside cells, this needs to be changed to a general boundary condition
            size_t const neighbor_j = neighbors[j];
            if(i < neighbor_j)
            {
                if(!tess.IsPointOutsideBox(neighbor_j))
                {
                    bool set_to_zero = false;
                    for(size_t k = 0; k < Nzero; ++k)
                        if(cells[neighbor_j].stickers[zero_indeces[k]])
                            set_to_zero = true;
                    double const Er_j = cells[neighbor_j].Erad * cells[neighbor_j].density * (set_to_zero ? zero_value : 1);
                    Vector3D const cm_ij = CM - tess.GetCellCM(neighbor_j);
                    Vector3D const grad_E = cm_ij * (1.0 / ScalarProd(cm_ij, cm_ij));
                    Vector3D const r_ij = point - tess.GetMeshPoint(neighbor_j);
                    double mid_D = 0.5 * (D[neighbor_j] + Dcell);
                    double const flux_limiter = flux_limiter_ ? CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j)) : 1;
                    mid_D *= flux_limiter;
                    double const flux = (self_zero || set_to_zero) ? tess.GetArea(faces[j]) * dt * CG::speed_of_light * 0.5 : ScalarProd(grad_E, r_ij * (tess.GetArea(faces[j]) / abs(r_ij))) * dt * mid_D; 
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
                    boundary_calc_.SetBoundaryValues(tess, i, neighbor_j, dt, cells, tess.GetArea(faces[j]), A[i][0], b[i], faces[j]);
            }
        }
    }
    for(size_t i = 0; i < Nlocal; ++i)
    {
        A[i].resize(max_neigh, 0);
        A_indeces[i].resize(max_neigh, max_size_t);
    }
}

void Diffusion::PostCG(Tessellation3D const& tess, std::vector<Conserved3D>& extensives, double const dt, std::vector<ComputationalCell3D>& cells,
        std::vector<double>const& CG_result)const
{
    size_t const N = tess.GetPointNo();
    bool const entropy = !(std::find(ComputationalCell3D::tracerNames.begin(), ComputationalCell3D::tracerNames.end(), std::string("Entropy")) ==
		ComputationalCell3D::tracerNames.end());
    size_t const entropy_index = static_cast<size_t>(std::find(ComputationalCell3D::tracerNames.begin(),
        ComputationalCell3D::tracerNames.end(), std::string("Entropy")) - ComputationalCell3D::tracerNames.begin());
    std::vector<size_t> zero_indeces;
    size_t const Nzero = zero_cells_.size();
    for(size_t i = 0; i < Nzero; ++i)
        zero_indeces.push_back(binary_index_find(ComputationalCell3D::stickerNames, zero_cells_[i]));
  
    for(size_t i = 0; i < N; ++i)
    {
        bool self_zero = false;
        for(size_t j = 0; j < Nzero; ++j)
            if(cells[i].stickers[zero_indeces[j]])
                self_zero = true;
        double const volume = tess.GetVolume(i);
        extensives[i].Erad = CG_result[i] * volume;
        double const T = cells[i].temperature;
        double const dE = fleck_factor[i] * CG::speed_of_light * dt * sigma_planck[i] * (CG_result[i] - T * T * T * T * CG::radiation_constant) * volume;
        extensives[i].energy += dE;
        extensives[i].internal_energy += dE;
        if(extensives[i].internal_energy < 0 || !std::isfinite(extensives[i].internal_energy))
            throw UniversalError("Bad internal energy in Diffusion::PostCG");
        cells[i].internal_energy = extensives[i].internal_energy / extensives[i].mass;
        cells[i].temperature = eos_.de2T(cells[i].density, cells[i].internal_energy, cells[i].tracers, ComputationalCell3D::tracerNames);
        cells[i].pressure = eos_.de2p(cells[i].density, cells[i].internal_energy, cells[i].tracers, ComputationalCell3D::tracerNames);
        cells[i].Erad = extensives[i].Erad / extensives[i].mass;
        if(entropy)
        {
            cells[i].tracers[entropy_index] = eos_.dp2s(cells[i].density, cells[i].pressure, cells[i].tracers, ComputationalCell3D::tracerNames);
            extensives[i].tracers[entropy_index] = cells[i].tracers[entropy_index] * extensives[i].mass;
        }

    }
}

void DiffusionSideBoundary::SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt, 
        std::vector<ComputationalCell3D> const& /*cells*/, double const Area, double& A, double &b, size_t const /*face_index*/)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        A += 0.5 * CG::speed_of_light * dt * Area;
        b += 2 * Area * dt * CG::stefan_boltzman * T_ * T_ * T_ * T_;
    }
}

void DiffusionSideBoundary::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
        E_outside = CG::radiation_constant * T_ * T_ * T_ * T_;
    else
        E_outside = new_E[index];
    v_outside = cells[index].velocity;
}

void DiffusionClosedBox::SetBoundaryValues(Tessellation3D const& /*tess*/, size_t const /*index*/, size_t const /*outside_point*/, double const /*dt*/, 
        std::vector<ComputationalCell3D> const& /*cells*/, double const /*Area*/, double& /*A*/, double& /*b*/, size_t const /*face_index*/)const
{}

void DiffusionClosedBox::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    E_outside = new_E[index];
    Vector3D normal = normalize(tess.GetMeshPoint(outside_point) - tess.GetMeshPoint(index));
    v_outside = cells[index].velocity;
    v_outside -= 2 * normal * ScalarProd(normal, v_outside);
}

double PowerLawOpacity::CalcDiffusionCoefficient(ComputationalCell3D const& cell) const
{
    return D0_ * std::pow(cell.density, alpha_) * std::pow(cell.temperature, beta_);
}

double PowerLawOpacity::CalcPlanckOpacity(ComputationalCell3D const& cell) const
{
    return CG::speed_of_light / (3 * CalcDiffusionCoefficient(cell));
}

void DiffusionXInflowBoundary::SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt, 
        std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        double const Er_j = left_state_.Erad * left_state_.density;
        double const Er = cells[index].Erad * cells[index].density;
        double const dx = tess.GetCellCM(index).x - tess.GetBoxCoordinates().first.x;
        Vector3D const grad_E = Vector3D(1, 0, 0) * (1.0 / (2 * dx));
        Vector3D const r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
        double mid_D = 0.5 * (D_calc_.CalcDiffusionCoefficient(cells[index]) + D_calc_.CalcDiffusionCoefficient(left_state_));
        double const flux_limiter = CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j));
        mid_D *= flux_limiter;
        double const flux = ScalarProd(grad_E, r_ij * (tess.GetArea(face_index) / abs(r_ij))) * dt * mid_D; 
        A += flux;
        b += flux * Er_j;
    }
    else
    {
        if(tess.GetMeshPoint(index).x < (tess.GetMeshPoint(outside_point).x - R * 1e-4))
        {
            double const Er_j = right_state_.Erad * right_state_.density;
            double const Er = cells[index].Erad * cells[index].density;
            double const dx = tess.GetBoxCoordinates().second.x - tess.GetCellCM(index).x;
            Vector3D const grad_E = Vector3D(-1, 0, 0) * (1.0 / (2 * dx));
            Vector3D const r_ij = tess.GetMeshPoint(index) - tess.GetMeshPoint(outside_point);
            double mid_D = 0.5 * (D_calc_.CalcDiffusionCoefficient(cells[index]) + D_calc_.CalcDiffusionCoefficient(right_state_));
            double const flux_limiter = CalcSingleFluxLimiter(grad_E * (Er - Er_j), mid_D, 0.5 * (Er + Er_j));
            mid_D *= flux_limiter;
            double const flux = ScalarProd(grad_E, r_ij * (tess.GetArea(face_index) / abs(r_ij))) * dt * mid_D; 
            A += flux;
            b += flux * Er_j;
        }
    }
}

void DiffusionXInflowBoundary::GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const
{
    double const R = tess.GetWidth(index);
    if(tess.GetMeshPoint(index).x > (tess.GetMeshPoint(outside_point).x + R * 1e-4))
    {
        E_outside = left_state_.Erad * left_state_.density;
        v_outside = left_state_.velocity;
    }
    else
    {
        if(tess.GetMeshPoint(index).x < (tess.GetMeshPoint(outside_point).x - R * 1e-4))
        {
            E_outside = right_state_.Erad * right_state_.density;
            v_outside = right_state_.velocity;
        }
        else
        {
            E_outside = new_E[index];
            Vector3D normal = normalize(tess.GetMeshPoint(outside_point) - tess.GetMeshPoint(index));
            v_outside = cells[index].velocity;
            v_outside -= 2 * normal * ScalarProd(normal, v_outside);
        }
    }
}