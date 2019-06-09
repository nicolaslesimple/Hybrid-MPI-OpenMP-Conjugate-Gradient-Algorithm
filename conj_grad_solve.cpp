#include "conj_grad_solve.hpp"
#include <omp.h>
#include <time.h>
#include <chrono>

/*****************************************************************************
                   Conjugate Gradient Algorithm cpp file

 Implementation of the main functions for the parallel Conjugate Gradient
 algorithm used in the main function.
*******************************************************************************/

using namespace std;
using vec = vector<double>;         // vector
using matrix = vector<vec>;         // matrix

void matrix_multiply_vector_openmp(const vector<vec> &sub_A, const vec &v, vec &result) {

    /***** Parallelization using OpenMP *****/
    /*
        In our code, when we have more that of proc meaning psize superior to 1, sub_A is a part of the matrix A (a subset of rows)
        defining Ax=b since we are 1D decomposing the matrix by rows

       #pragma omp simd : Applied to a loop to indicate that the loop can be transformed into a SIMD loop.
        This is the case.
       #pragma omp for simd : Specifies that a loop that can be executed concurrently using SIMD instructions,
        and that those iterations will also be executed in parallel by threads in the team.
       #pragma omp parallel for simd : Shortcut for specifying a parallel construct containing one
        for simd construct and no other statements.
    */

    size_t sub_size = sub_A.size();
    //#pragma omp simd // This command tells use that we can use simd loop.
    #pragma omp parallel for simd // This command allows to divide the work.
    for (size_t i = 0; i < sub_size; i++)
        result[i] = dot_product(sub_A[i], v);  // Here we call another parallelized function for the dot product
}

void vec_lin_combo(double a, const vec &u, double b, const vec &v, vec &result) {

    /***** Without Parallelization *****/
    /*
        In this function, none of the operation is parallelized.
        In fact, it seems after experimenting our code that it is faster without this pragma call
        This true at least for the matrices size we are currently working on.
        This observation can be explained by the fact that overhead of creating threads is larger than operation.
    */

    size_t n = u.size();
    for (size_t j = 0; j < n; j++)
        result[j] = a * u[j] + b * v[j];
}

void vec_lin_combo_openmp(double a, const vec &u, double b, const vec &v, vec &result) {

     /***** Parallelization using OpenMP *****/
     /*
        In this function, the linear combination operation is parallelized with OpenMP.

       #pragma omp parallel for : Shortcut for specifying a parallel construct containing one
        or more associated loops and no other statements.
    */

    size_t n = u.size();
    #pragma omp parallel for // This command allows to divide the work.
    for (size_t j = 0; j < n; j++)
        result[j] = a * u[j] + b * v[j];
}

double dot_product_normal(const vec &u, const vec &v){

    // None parallelization done in this function.
    size_t n = u.size();
    double value = 0;
    for (size_t i = 0; i < n; i++) {
        value += u[i] * v[i];
    }
    return value;
}

double dot_product(const vec &u, const vec &v) {

    return __gnu_parallel::inner_product(u.begin(), u.end(), v.begin(), 0.0);  // parallelized inner product
}


double dot_product_mpi(const vec &sub_u, const vec &sub_v) {

     /***** Parallelization using OpenMP and MPI*****/
     /*
        In this function, the dot product operation is parallelized with OpenMP.
        It performs a reduction over the sub-vectors which are passed to it.
        In fact, All_Reduce broadcasts the value to all psize.

       #pragma omp parallel for : Shortcut for specifying a parallel construct containing one
        or more associated loops and no other statements.
       #reduction(reduction-identifier:list) : Specifies a reduction-identifier and one or more list items.
        The reduction-identifier must match a previously declared reduction-identifier of the same name and
        type for each of the list items.
       # MPI_Allreduce : Combines values from all processes and distributes the result back to all processes.
        This routine is threads safe.
    */

    // Initialization
    double product;
    size_t length = sub_u.size();
    double sub_prod = 0.0;

    // Parallelization with OpenMP and MPI
    #pragma omp parallel for reduction(+:sub_prod)
    for (size_t i = 0; i < length; i++) {
        sub_prod += sub_u[i] * sub_v[i];
    }
    MPI_Allreduce(&sub_prod, &product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // We do a reduction over sub_prod to get the total dot product
    return product;
}


double vector_norm(const vec &v) {
    // We use the dot_product_mpi function and we apply on it the square root to have the norm
    return sqrt(dot_product(v, v));
}


vec conj_grad_solver(const matrix &sub_A, const vec &b, const double tolerance, const vec &initial_guess, int &total_iters) {

    /***** Parallelization using OpenMP and MPI *****/
    /*
        Our goal in this algorithm is to find the x vector knowing matrix A and vector b such that Ax=b.
        In our code, when we have more that of proc meaning psize superior to 1, sub_A is a part of the
        matrix A (a subset of rows) defining Ax=b.

        In addition to matrix A and vector b, our function take as argument the tolerance which is a double
        defining when the algorithm succeed to converge.
        Moreover, it takes the initial guess corresponding to the initial solution.

        In this function we will mainly use MPI command line. However, we will call functions that are using
        OpenMP command.

        #MPI_Allgatherv : Gathers data from all tasks and deliver the combined data to all tasks.
        #MPI_Gatherv : Gathers into specified locations from all processes in a group.

        Just one last thing to be careful : make sure the system size is big enough for the number of processors.
    */

    /********* Initialisation of MPI **********/
    int psize, prank;
    MPI_Comm_size (MPI_COMM_WORLD, &psize);
    MPI_Comm_rank (MPI_COMM_WORLD, &prank);

    /********* Decomposition and Initialisation of the variables needed to CG algo **********/
    int row_cnt[psize];  // This variable allows us to have the psize for each process (i.e number of rows in each prank).
                          // In fact, when the size of the system is not evenly divisible, this variable can be different.
    int row_disp[psize]; // This variable will be used in one of the final MPI command which is gatherv.
                          // It corresponds to the displacement from start of vector.
    size_t m = b.size();

    for (int i = 0; i < psize; i++) {
        row_cnt[i] = m/psize;
        row_disp[i] = i * m/psize;  // In our algorithm, prank 0 is also doing work. Otherwise we should have done (psize - 1).
    }
    // Initialization of needed variables
    int sub_vector_size = m/static_cast<size_t>(psize);
    vec sub_x(sub_vector_size); // After Gathering, this vector will be our answer to the problem.
    vec sub_r(sub_vector_size); // R is the residual and it will say if we succeed to converge.
    // Here we want to create the corresponding subpart of of r and x to sub_A given as input.
    // In fact, in the main the matrix A was already divided. We do the same thing for r and x.
    for (size_t i = 0; i < sub_vector_size; i++) {
        int position = (sub_vector_size)*static_cast<size_t>(prank) + i;
        sub_r[i] = b[position];
        sub_x[i] = initial_guess[position];
    }

    /********** Implementation of the algorithmic part of the Conjugate Gradient Method **********/
    /*              To see more details on the algorithmic part, look at the report              */
    /*********************************************************************************************/

    /* Initialization */
    vec x(m);  // This vector will be the solution of our problem. It will be return to the main.
    vec p = b;  // Here we did not decomposed the vector for now.
    vec sub_p = sub_r; // Here we decomposed the previous vector according to prank and psize. P is the search direction.
    vec sub_r_old;
    vec result1(sub_vector_size), result2(sub_vector_size), result3(sub_vector_size);
    vec vec_sub_A_multiply_by_p(sub_vector_size), vec_sub_A_multiply_by_sub_x(sub_vector_size);
    int max_iter = 100000;
    double alpha_num, denominator_alpha, alpha, alpha_num_old, beta;
    alpha_num = dot_product_mpi(sub_r, sub_r); // We do this first calculation before the for loop as it is just for the first iteration

    /* Main part of CG algo */
    for (int i = 0; i < max_iter; i++) {

        // Calculate all variable to find alpha
        sub_r_old = sub_r;
        alpha_num_old = alpha_num;
        matrix_multiply_vector_openmp(sub_A, p, vec_sub_A_multiply_by_p);
        denominator_alpha = dot_product_mpi(sub_p, vec_sub_A_multiply_by_p);
        // Alpha calculation
        alpha = alpha_num/denominator_alpha;
        // Update of the solution x and the residual r
        vec_lin_combo(1.0, sub_x, alpha, sub_p, result1);
        sub_x = result1;
        vec_lin_combo(1.0, sub_r, -alpha, vec_sub_A_multiply_by_p, result2);
        sub_r = result2;
        alpha_num = dot_product_mpi(sub_r, sub_r);

        // Here we test the convergence of our algorithm
        if (vector_norm(sub_r) < tolerance) {//sqrt(alpha_num)
            if (prank == 0) {
                cout << "Converged at iter = " << i << endl; // Print the convergence time on user interface.
            }
            total_iters = i;
            break; // We need to be careful with a break in a parallelization algorithm. Here it's working well as all threads are merged.
        }

        // Beta parameter update
        beta = alpha_num/alpha_num_old;
        // Update the search direction p
        vec_lin_combo(1.0, sub_r, beta, sub_p, result3);
        sub_p = result3;

        // We need to take back the entire vector p as we need it in the matrix_multiply_vector_openmp for the next iteration. Thus we gather the vector.
        MPI_Allgatherv(&sub_p.front(), row_cnt[prank], MPI_DOUBLE, &p.front(), row_cnt, row_disp, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    /* Return the solution to the main after convergence check */
    // Again we gather all processes to obtain the final vector x that will be return to the main function.
    MPI_Gatherv(&sub_x.front(), row_cnt[prank], MPI_DOUBLE, &x.front(), row_cnt, row_disp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return x;
}
