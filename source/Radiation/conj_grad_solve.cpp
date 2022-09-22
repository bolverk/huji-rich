#include "conj_grad_solve.hpp"
#include <limits>

namespace CG
{
    // Matrix times vector
    void mat_times_vec(const mat &sub_A_values, const size_t_mat &sub_A_indices, const std::vector<double> &v, 
        std::vector<double> &result)
    {

        // NOTE: when using MPI with > 1 proc, A will be only a sub-matrix (a subset of rows) of the full matrix
        // since we are 1D decomposing the matrix by rows

        size_t const sub_num_rows = sub_A_values.size();
        size_t const sub_num_cols = sub_A_values[0].size();
        result.resize(sub_num_rows, 0);
        double dot_prod;
        for (size_t i = 0; i < sub_num_rows; i++) {
            dot_prod = 0;  // rezero the dot_prod buffer. we need this buffer so we can make it private to the thread to avoid race conditions.
            for (size_t j = 0; j < sub_num_cols; j++) {
                if(sub_A_indices[i][j] == max_size_t)
                    break;
                dot_prod += sub_A_values[i][j] * v[sub_A_indices[i][j]]; 
            }
            result[i] = dot_prod;
        }
    }

    std::vector<double> vector_rescale(std::vector<double> const& a, std::vector<double> const& b)
    {
        size_t const N = a.size();
        if(a.size() != b.size())
            throw UniversalError("Sizes do not match in vector_rescale");
        std::vector<double> res(N);
        for(size_t i = 0; i < N; ++i)
            res[i] = a[i] * b[i];
        return res;
    }

    // Linear combination of vectors
    void vec_lin_combo(double a, const std::vector<double> &u, double b, const std::vector<double> &v, 
        std::vector<double> &result)
    {
        if(u.size() != v.size())
            throw UniversalError("Unequal vector sizes in vec_lin_combo");
        size_t n = u.size();
        result.resize(n, 0);
        for (size_t j = 0; j < n; j++)
            result[j] = a * u[j] + b * v[j];
    }


    // performs a reduction over the sub-vectors which are passed to it... All_Reduce broadcasts the value to all procs
    double mpi_dot_product(const std::vector<double> &sub_u, const std::vector<double> &sub_v) // need to pass it the buffer where to keep the result
    {
        if(sub_u.size() != sub_v.size())
            throw UniversalError("Unequal vector sizes in vec_lin_combo");
        size_t length = sub_u.size();

        double sub_prod = 0.0;
        for (size_t i = 0; i < length; i++) {
            sub_prod += sub_u[i] * sub_v[i];
        }
#ifdef RICH_MPI
        // do a reduction over sub_prod to get the total dot product
        MPI_Allreduce(MPI_IN_PLACE, &sub_prod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        return sub_prod;
    }
    
    void build_M(const mat &sub_A_values, const size_t_mat &sub_A_indices, std::vector<double> &M)
    {
        size_t const sub_num_rows = sub_A_values.size();
        size_t const sub_num_cols = sub_A_values[0].size();
        M.resize(sub_num_rows);
        for (size_t i = 0; i < sub_num_rows; i++) {
            for (size_t j = 0; j < sub_num_cols; j++) {
                if(sub_A_indices[i][j] == max_size_t)
                    break;
                if(sub_A_indices[i][j] == i)
                {
                    M[i] = 1.0 / sub_A_values[i][j];
                    //M[i] = 1.0;
                    break;
                }
            }
        }
    }

    static bool abs_compare(double a, double b)
    {
        return (std::abs(a) < std::abs(b));
    }

    std::vector<double> conj_grad_solver(const double tolerance, int &total_iters,
        Tessellation3D const& tess, std::vector<ComputationalCell3D> const& cells,
        double const dt, MatrixBuilder const& matrix_builder, double const time)  //total_iters is to store # of iters in it
    {
        size_t Nlocal = tess.GetPointNo();
        
        // NOTE: when using MPI with > 1 proc, A will be only a sub-matrix (a subset of rows) of the full matrix
        // since we are 1D decomposing the matrix by rows
        // b will be the full vector

        int nprocs = 1, rank = 0;
    #ifdef RICH_MPI
        MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    #endif
    
        int max_iter = 1000;

        mat A;
        size_t_mat A_indeces;
        std::vector<double> b;
        std::vector<double> sub_x; // this is for the initial guess
        matrix_builder.BuildMatrix(tess, A, A_indeces, cells, dt, b, sub_x, time);
        std::vector<double> M; // The preconditioner
        build_M(A, A_indeces, M);

        std::vector<double> r_old, sub_a_times_p;
        std::vector<double> sub_r;
#ifdef RICH_MPI
        MPI_exchange_data2(tess, sub_x, true);
#endif
        mat_times_vec(A, A_indeces, sub_x, sub_a_times_p);
        // Find maximum value of A, this is used for normalization of the error
        double maxA[2] = {0, 0};
        for(size_t i = 0; i < A.size(); ++i)
        {
            maxA[0] = std::max(maxA[0], std::abs(sub_x[i]));
            maxA[1] = std::max(maxA[1], std::abs(b[i]));
        }       
#ifdef RICH_MPI
        MPI_Allreduce(MPI_IN_PLACE, &maxA, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
        vec_lin_combo(1.0, b, -1.0, sub_a_times_p, sub_r);    
        std::vector<double> sub_p(sub_r);
        sub_p.resize(Nlocal);
        sub_x.resize(Nlocal);
        sub_p = vector_rescale(sub_p, M);
        std::vector<double> result1, result2, result3, p(sub_p);
        std::vector<double> old_x = sub_x;
        size_t Ntotal = Nlocal;
#ifdef RICH_MPI
        MPI_exchange_data2(tess, p, true);
        MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
        double sub_r_sqrd = mpi_dot_product(sub_r, sub_p);
        double sub_r_sqrd_convergence = mpi_dot_product(sub_r, sub_r);
        double sub_r_sqrd_old = 0, sub_p_by_ap = 0, alpha = 0, beta = 0;
        bool good_end = false;
        double max_values[2] = {0, 0};
        // Main Conjugate Gradient loop
        // this loop must be serial b/c CG is an iterative method
        for (int i = 0; i < max_iter; i++) {
            // note: make sure matrix is big enough for the number of processors you are using!
            r_old = sub_r;                 // Store previous residual
            sub_r_sqrd_old = sub_r_sqrd;  // save a recalculation of r_old^2 later

            mat_times_vec(A, A_indeces, p, sub_a_times_p);  //split up with MPI and then finer parallelize with openmp

            sub_p_by_ap = mpi_dot_product(sub_p, sub_a_times_p);

            alpha = sub_r_sqrd / (sub_p_by_ap + std::numeric_limits<double>::min() * 100);         

            // Next estimate of solution
            vec_lin_combo(1.0, sub_x, alpha, sub_p, result1);
            sub_x = result1;
            vec_lin_combo(1.0, sub_r, -alpha, sub_a_times_p, result2);
            sub_r = result2;
            size_t max_loc = 0;
            max_values[1] = 0;
            max_values[0] = 0;
            for(size_t j = 0; j < Nlocal; ++j)
            {
                double const local_scale = std::abs(b[j]);
                if(std::abs(sub_r[j]) / (local_scale + maxA[1] * 1e-3) > max_values[1])
                    max_loc = j;
                max_values[1] = std::max(max_values[1], std::abs(sub_r[j]) / (local_scale + maxA[1] * 0.1));
                max_values[0] = std::max(max_values[0], std::abs(sub_x[j] - old_x[j]) / (std::abs(sub_x[j]) + std::numeric_limits<double>::min() * 100 + maxA[0] * 1e-4));
            }
            sub_r_sqrd_convergence = mpi_dot_product(sub_r, sub_r);
            result2 = vector_rescale(sub_r, M);
            sub_r_sqrd = mpi_dot_product(sub_r, result2);

    #ifdef RICH_MPI
            MPI_Allreduce(MPI_IN_PLACE, max_values, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    #endif
            // recall that we can't have a 'break' within an openmp parallel region, so end it here then all threads are merged, and the convergence is checked
            // Convergence test
            if (std::sqrt(sub_r_sqrd_convergence / Ntotal) < tolerance * maxA [1]
                && max_values[1] < 1e-6 && max_values[0] < 1e-6) { // norm is just sqrt(dot product so don't need to use a separate norm fnc) // vector norm needs to use a all reduce!
#ifdef RICH_MPI
                if(rank == 0)
#endif
                std:: cout << "Converged at iter = " << i << std::endl;
                total_iters = i;
                good_end = true;
                break;
            }
            old_x = sub_x;
            beta = sub_r_sqrd / sub_r_sqrd_old;       
            
            vec_lin_combo(1.0, result2, beta, sub_p, result3);             // Next gradient
            sub_p = result3;
            p = sub_p;
    #ifdef RICH_MPI
            MPI_exchange_data2(tess, p, true);
    #endif
        }
        if(not good_end)
            throw UniversalError("CG did not converge");
#ifdef RICH_MPI
        MPI_exchange_data2(tess, sub_x, true);
#endif
        return sub_x;
    }

}