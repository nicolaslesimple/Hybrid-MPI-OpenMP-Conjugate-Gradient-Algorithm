# Definition of the compiler option, libraries and the target
FLAGS = -O3 -qopt-report -Wall -fopenmp -pg -Wextra

# Creation of the executable file
SRCS =   main.cpp conj_grad_solve.cpp
OBJS = $(subst .cpp,.o,$(SRCS))
all: cg

# Creation of the object file -o
cg: $(OBJS)  
	mpiicc ${FLAGS} -o cg $(OBJS)
	export OMP_NUM_THREADS=1  # use this for OpenMP

%.o: %.cpp 
	mpiicc $(FLAGS) -c $<

# Clean
clean:
	rm *.o ./cg
