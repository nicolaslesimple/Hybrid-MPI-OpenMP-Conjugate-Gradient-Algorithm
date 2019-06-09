#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <parallel/numeric> // for gnu_inner_produc

/*****************************************************************************
                   Conjugate Gradient Algorithm hpp file

 Implementation of the main functions for the parallel Conjugate Gradient
 algorithm used in the main function.
*******************************************************************************/

using namespace std;
using vec = vector<double>;         // vector
using matrix = vector<vec>;         // matrix

/* Matrix times vector product */
void matrix_multiply_vector_openmp(const vector<vec> &A, const vec &v, vec &result);

/* Linear combination Operation without Parallelization */
void vec_lin_combo(double a, const vec &u, double b, const vec &v, vec &result);

/* Linear combination Operation with OpenMP */
void vec_lin_combo_openmp(double a, const vec &u, double b, const vec &v, vec &result);

/* Dot product without parallelization  */
double dot_product(const vec &u, const vec &v);

/* Dot product of vectors with OpenMP and MPI */
double dot_product_mpi(const vec &sub_u, const vec &sub_v);

/* Norm of a vector */
double vector_norm(const vec &v);

/* Conjugate Gradient Algorithm Implementation */
vec conj_grad_solver(const matrix &A, const vec &b, const double tolerance, const vec &initial_guess, int &total_iters);
