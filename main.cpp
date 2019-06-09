#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <mpi.h>
#include <time.h>
#include <chrono>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <iomanip>
#include "conj_grad_solve.hpp"

/****************************************************************************************************
                             Conjugate Gradient Algorithm main file

 This file allows the usage of our parallel implementation of the conjugate gradient algorithm using
 a hybrid of  distributed (MPI) and shared (OpenMP) memory approach for dense matrices.
***************************************************************************************************/

using namespace std;
using vec = vector<double>;         // Vector structure we will use in the entire code
using matrix = vector<vec>;         // Matrix structure we will use in the entire code


/* Declaration of the function used to print a vector on the user interface */
void print(const vec &x) {
   size_t N = x.size();
   for (size_t i = 0; i < N; i++) {
      cout << fixed << setprecision(10) << x[i] << '\n';
   }
   cout<< '\n';
}

/* Declaration of the function used to print a matrix on the user interface */
void print(const matrix &A) {
   size_t M = A.size();
   size_t N = A[0].size();
   for (size_t i = 0; i < M; i++) {
      for (size_t j = 0; j < N; j++) {
         cout << fixed << setw(10) << setprecision(5) << A[i][j];
      }
      cout << '\n';
   }
}


int main(int argc, char **argv) {


    /********** Initialisation **********/
    double tolerance = 1e-12;
    double tolerance_error = 1e-5;
    int num_solves = 5;  // This integer represents the number of solves of CG to do. It will give better statistics of cpu time.
    int psize, prank;
    int total_iters;

    // Matrix size
    int M, N;
	M = atoi(argv[1]); // Take the size of the system with the input given by the user in the command line
	N = M; // We are working with square matrix

    // We initialize the several elements defining the system Ax=b
    matrix A(N, vector<double> (N));
    vec b(N,0), x;
    // We initialize the several vectors needed for the CG algorithm of the size of the system
    vec initial_guess(N,0), r(N,0), p(N,0);

    // MPI Initialisation
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &psize);
    MPI_Comm_rank (MPI_COMM_WORLD, &prank);

    /*********** Here we read the previously created matrix and vector
    and create the corresponding variables defining our system ***********/
    ifstream matrixfile("matrix");
    if(!(matrixfile.is_open())) {
        cout<<"Error: file not found"<<endl; // Check if the file is really open
        return 0;
    }

    // Attribute value to the file to matrix A
    for(int i = 0; i < M; i++){
        vec tmp;
        for(int j = 0; j < N; j++){
            matrixfile >> A[i][j];
        }
    }
    // Attribute value to the vector b
    for(int i = 0; i < M; i++){
        matrixfile >>b[i];
        initial_guess[i] = 0;
    }
    matrixfile.close();

    /*********** Creation of different subparts of the matrix A according to the
    psize defined for the simulation ***********/
    matrix sub_A(N/static_cast<size_t>(psize), vector<double> (N));
    for (size_t i = 0; i < N/static_cast<size_t>(psize); i++) {
        for (size_t j = 0; j < N; j++) {
            sub_A[i][j] = A[static_cast<size_t>(prank) * N/static_cast<size_t>(psize) + i][j];
        }
    }

    /***********  Repeat Conjugate gradient algorithm a few times. This process will
    allow to obtain more reliable statistics of CPU time. ***********/
    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now(); // Measure of the time taken for the entire for loop
    for (int i = 0; i < num_solves; i++) {
        x = conj_grad_solver(sub_A, b, tolerance, initial_guess, total_iters); // We store the results on the previously defined
    }


    /*********** Measure of the time and displaying the results on the user interface ***********/
    chrono::high_resolution_clock::time_point finish = chrono::high_resolution_clock::now();
    chrono::duration<double> time = chrono::duration_cast<chrono::duration<double>>(finish-start);

    if (prank == 0) {
        // We display the system Ax=b --> decomment the code to allow the printing
        /*cout << "Matrix A: " << endl;
        print(A);
        cout << "Vector b: " << endl;
        print(b);
        // We display our found solution
        cout << "Solution x: " << endl;
        print(x);
        // We check if our solution make sense by multiplying it with A as Ax = b
        cout << "Check A*x = b :" << endl;*/
        vec A_times_x(x.size());
        matrix_multiply_vector_openmp(A, x, A_times_x);
        //print(A_times_x);

        // We will compute the error between the solution we found multiply by A and between b vector
        double mean_error = 0;
        vec error(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            error[i] = abs(A_times_x[i] - b[i]); // The error vector is fill up with the difference for each points
            mean_error += error[i];
        }
        if (*max_element(error.begin(), error.end()) > tolerance_error)
            cout << "Error in solution is larger than " << tolerance_error << endl;
            cout << "In fact, the error is :" << mean_error/N <<endl;

        // Calculation of the different time taken by our simulation
        double total_cpu_time = time.count();
        double cpu_time = total_cpu_time/num_solves;
        double cpu_time_per_iter = cpu_time/total_iters;
        // Print all these values on the user interface
        cout << " Total CPU time = " << total_cpu_time << endl;
        cout << " CPU time per CG solve = " << cpu_time << endl;
        cout << " CPU time per iter = " << cpu_time_per_iter << endl;
    }

    // Finalize MPI processes
    MPI_Finalize();
}
