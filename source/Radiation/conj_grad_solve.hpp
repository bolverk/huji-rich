#ifndef CG_HPP
#define CG_HPP 1

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#ifdef RICH_MPI
#include <mpi.h>
#include "source/mpi/mpi_commands.hpp"
#endif
#include <boost/container/small_vector.hpp>
#include "source/3D/GeometryCommon/Tessellation3D.hpp"
#include "source/newtonian/three_dimensional/computational_cell.hpp"
#include "source/newtonian/three_dimensional/conserved_3d.hpp"
#include "source/misc/utils.hpp"

namespace CG
{
    using vec = boost::container::small_vector<double, 20>;         // vector
    using vec_size_t = boost::container::small_vector<size_t, 20>;         // vector
    using mat = std::vector<vec>;            // matrix (=collection of (row) vectors)
    using size_t_mat = std::vector<vec_size_t>;
    double constexpr speed_of_light = 2.99792458e10;
    size_t constexpr max_size_t = std::numeric_limits<size_t>::max();
    double constexpr stefan_boltzman = 5.670374e-5;
    double constexpr radiation_constant = 4 * stefan_boltzman / speed_of_light;

    //! \brief Class that build the data for the solution of the linear system A*x=b
    class MatrixBuilder
    {
         public:
            MatrixBuilder(std::vector<std::string> const zero_cells = std::vector<std::string> ()) : zero_cells_(zero_cells){}
        /*!
            \brief Builds the initial conditions for the CG method to solve A*x=b
            \param tess The tessellation
            \param A The A matrix to build
            \param A_indeces The indeces of the values in A, this is needed since A is sparse
            \param cells The computational cells
            \param dt The time step
            \param b The b vector to calculate
            \param x0 The initial solution guess
            \param current_time The time
        */
        virtual void BuildMatrix(Tessellation3D const& tess, mat& A, size_t_mat& A_indeces, std::vector<ComputationalCell3D> const& cells,
            double const dt, std::vector<double>& b, std::vector<double>& x0, double const current_time) const = 0;
        /*!
        \brief This method does post processing after the CG has finished (e.g. update the thermal energy)
        \param tess The tesselation
        \param extensives The extensives
        \param dt The time step
        \param cells The primitive variables
        \param CG_result The result from the CG
        */
        virtual void PostCG(Tessellation3D const& tess, std::vector<Conserved3D>& extensives, double const dt, std::vector<ComputationalCell3D>& cells,
            std::vector<double>const& CG_result)const = 0;

        std::vector<std::string> const zero_cells_;
    };

    //! The fastest implementation of conjugate gradient algorithm, using data-based parallelism only
    std::vector<double> conj_grad_solver(const double tolerance, int &total_iters,
        Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells,
        double const dt, MatrixBuilder const& matrix_builder, double const time);
}

#endif