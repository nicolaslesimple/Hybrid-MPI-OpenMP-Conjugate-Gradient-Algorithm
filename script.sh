#!/bin/bash

module load intel
module load intel-mpi

make

python Random_matrix.py 1000
srun -t1000  -n1 ./cg 1000


