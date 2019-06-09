import random
import sys
import numpy as np

"""
This function allows the generation of a random dense square matrix and two vectors.
It creates the matrix A, the vector b and x corresponding to the linear equation Ax = b.
This matrix will be used to test our parallelized Conjugate Gradient Algorithm.
To generate the vector, the user need to use the following command line in the terminal where N = 10 :
>>> python Random_matrix.py 10
"""

if (len(sys.argv) > 1):

    # Take and save the desired size of rhe matrix given by the user
    M = int(sys.argv[1])
    N = M

    # Initialisation of the matrix and the two vector
    x = [0 for i in range(M)] # x correspond to the x vector we should find with CG algorithm
    b = [0 for i in range(M)] # b is the vector in the following equation Ax = b
    A = [[0 for i in range(N)] for j in range(M)] # A is the matrix defining the following system Ax = b

    # Random initialisation of the vector x
    for i in range(M):
        x[i] = random.random()

    # These for loops allows creation of matrix symmetric definite positive A
    for i in range(M):
        for j in range(N):
            if i <= j:
                A[i][j] = random.random() # Gives random numbers
                if i < j:
                    if random.random() < 0.1: # Allows to have 0 in the matrix
                        A[i][j] = 0
                A[j][i] = A[i][j] # Allows the symmetrical property
                b[i] += (A[i][j] * x[j]) # Give the right number to vector b to allows a right definition of Ax = b

    A = np.dot(A,A)
    b = np.dot(A,b)

    # Saving part of the function for A
    f = open('matrix','w') # Open the file where the matrix will be saved
    for i in range(M):
        for j in range(N):
            f.write(str(A[i][j])) # Write value of the matrix inside the file
            f.write('\t')
        f.write('\n')

    # Saving of vector b in the same file
    for i in range(M):
        f.write(str(b[i]))
        f.write('\n')
    f.close()

    # Saving of vector x in a file called x
    f = open('x', 'w')
    for i in range(M):
        f.write(str(x[i]))
        f.write('\n')
    f.close()

else:
    # If the size of the wanted system was not given by the user, a message is print on the interface.
    print ("please input a number")

