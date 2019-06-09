# Instruction to run the Paralelized implementation of the Conjugate Gradient Descent on EPFL server :

* Connect on deneb1 server:
    
    `ssh identifier@deneb1.epfl.ch`
    
***Then, the user have two choices : do it manually or use the bash script.sh.***

#### With the Bash Script :

* Just use the following command on the terminal inside the folder where the bash script is :

    `./script.sh`
    
If the user wants to choose the size of the system of linear equation he wants to solve, he just have to open the bash script and change the N value indicated below.
In addition, if he wants to change the number of distributed procs for MPI, he have to change the number p in the bash script.

    python Random_matrix.py N
    
    srun  -n p ./cg N
    
#### Do it manually :

* First, these particular modules are needed in order to be able to use :

    `module purge`
    
    `module load intel`
    
    `module load intel-mpi`
    
* Then to compile the code the following command line need to be used :

    `make`

    NOTE: the number of openmp threads is controlled by export OMP_NUM_THREADS= which is set inside of the makefile.
    
    We can also change the num threads environmental variable intereactively by putting into the terminal: export OMP_NUM_THREADS=[desired number of threads]

* In a third time, the system of equation we want to solve with our Conjugate Gradient algortihm need to be created.
To do that the user need to choose the size of the system. Then one random symetric definite positive matrix and two vector will be created to defined a system of the following form : A*x=b.
In the command line, the N represent the size that can be choosen by the user.

    `python Random_matrix.py N`
    
* Finally, we can run our parallelized implementation of the Conjugate Gradient Method with the following command :
To work well, the N integer should correspond to the one defined in the previous command.
In addition, the number p can also be choosen by the user. It represents the number of distributed procs for MPI.   
    
    `srun  -n p ./cg N`

* IMPORTANT NOTE : To use our implementation, the user need to verify the number of processes with the sizeof the system n.  In fact, to allow a well behavior of our code, we assume that the number of elements n are greater than or equal to number of processes p.


    

