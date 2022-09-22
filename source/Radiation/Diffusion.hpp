#ifndef DIFFUSION_HPP
#define DIFFUSION_HPP 1

#include "conj_grad_solve.hpp"
#include "source/newtonian/common/equation_of_state.hpp"

using namespace CG;

namespace CG
{
     double CalcSingleFluxLimiter(Vector3D const& grad, double const D, double const cell_value);
}

//! \brief Abstract class for calculating the needed data for diffusion
class DiffusionCoefficientCalculator
{
public:
/*!
    \brief Calculates the diffusion coefficient
    \param cell The primitive variables
    \return The diffusion coefficient (default units are cm^2/sec)
*/
virtual double CalcDiffusionCoefficient(ComputationalCell3D const& cell) const = 0;
/*!
    \brief Calculates the Planck opacity
    \param cell The primitive variables
    \return The planck opacity (default units are 1/cm)
*/
virtual double CalcPlanckOpacity(ComputationalCell3D const& cell) const = 0;
};

//! \brief Class for assigning boundary conditions for diffusion
class DiffusionBoundaryCalculator
{
    public:
    /*!
\brief Sets the boundary values for the matrix build in A*x=b
\param tess The tesselation
\param index The index of the cell that is adjacent to a boundary
\param outside_point The index of the cell that is outside
\param cell The primitve cells
\param A the value in the A matrix to change, given as input and output
\param b the value in the b vector to change, given as input and output
\param Area The area of the interface between the two cells
\param dt The time step
\param face_index The index of the interface
    */
virtual void SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
    std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const = 0;
    /*!
\brief Sets the outside values for the Froce calcualtion if needed
\param tess The tesselation
\param index The index of the cell that is adjacent to a boundary
\param outside_point The index of the cell that is outside
\param cell The primitve cells
\param E_outside The value of the energy in the outside cell, given as output
\param v_outside The value of the velocity, given as output
\param new_E The values of the new energy after the CG step
    */
virtual void GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
    std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const = 0;
};

//! \brief Class with constant blackbody temperature on the left x side and zero flux on other sides
class DiffusionSideBoundary : public DiffusionBoundaryCalculator
{
    public:
    /*!
    \brief Class constructor
    \param T Boundary temperature, in kelvin
    */
    DiffusionSideBoundary(double const T): T_(T){}

    void SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const override;

    void GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
        std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const override;
    private:
        double const T_;
};

//! \brief Class with constant blackbody temperature on the left x side and zero flux on other sides
class DiffusionClosedBox : public DiffusionBoundaryCalculator
{
    public:
    void SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const override;
    
    void GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
        std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const override;
};

//! \brief Class for calculating diffusion matrix data for the CG solver
class Diffusion : public CG::MatrixBuilder
{
public:
/*!
\brief Class constructor
\param D_coefficient_calc Class for calcualting the diffusion coefficients
\param eos The equation of state
\param boundary_calc Class to calcualte the values for the boundary conditions
*/
    Diffusion(DiffusionCoefficientCalculator const& D_coefficient_calc, DiffusionBoundaryCalculator const& boundary_calc,
        EquationOfState const& eos, std::vector<std::string> const zero_cells = std::vector<std::string> (), bool const flux_limiter = true) : D_coefficient_calcualtor(D_coefficient_calc),
        boundary_calc_(boundary_calc), eos_(eos), flux_limiter_(flux_limiter), sigma_planck(), fleck_factor(), CG::MatrixBuilder(zero_cells) {}

    void BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells, 
            double const dt, std::vector<double>& b, std::vector<double>& x0, double const current_time) const override;

    void PostCG(Tessellation3D const& tess, std::vector<Conserved3D>& extensives, double const dt, std::vector<ComputationalCell3D>& cells,
        std::vector<double>const& CG_result)const override;

    DiffusionCoefficientCalculator const& D_coefficient_calcualtor;
    mutable std::vector<double> sigma_planck, fleck_factor, D;
    DiffusionBoundaryCalculator const& boundary_calc_;
    bool const flux_limiter_;
private:  
    EquationOfState const& eos_;   
};

//! D=D0*rho^alpha*T^beta
class PowerLawOpacity: public DiffusionCoefficientCalculator
{
private:
    double const D0_, alpha_, beta_;
public:
    PowerLawOpacity(double const D0, double const alpha, double const beta): D0_(D0), alpha_(alpha), beta_(beta){}

    double CalcDiffusionCoefficient(ComputationalCell3D const& cell) const override;

    double CalcPlanckOpacity(ComputationalCell3D const& cell) const override;
};
//! \brief Class with constant states on the x sides and zero flux on other sides
class DiffusionXInflowBoundary : public DiffusionBoundaryCalculator
{
    public:
    /*!
    \brief Class constructor
    \param T Boundary temperature, in kelvin
    */
    DiffusionXInflowBoundary(ComputationalCell3D const& left_state, ComputationalCell3D const& right_state,
        DiffusionCoefficientCalculator const& D_calc): left_state_(left_state), right_state_(right_state), D_calc_(D_calc){}

    void SetBoundaryValues(Tessellation3D const& tess, size_t const index, size_t const outside_point, double const dt,
        std::vector<ComputationalCell3D> const& cells, double const Area, double& A, double &b, size_t const face_index)const override;
    
    void GetOutSideValues(Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells, size_t const index, size_t const outside_point,
        std::vector<double> const& new_E, double& E_outside, Vector3D& v_outside)const override;
    private:
        ComputationalCell3D const& left_state_, right_state_;
        DiffusionCoefficientCalculator const& D_calc_;
};
#endif