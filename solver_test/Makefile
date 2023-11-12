FC=gfortran

# Flags for g95 that you might also want to use with gfortran
FFLAGS=-g -pedantic -Wall -fbounds-check 

all: test_solver
	./test_solver

solver_module.o: solver_module.f90
	$(FC) -o solver_module.o -c solver_module.f90

test_solver.o: test_solver.f90 solver_module.o
	$(FC) $(FFLAGS) -c test_solver.f90 -o test_solver.o

test_solver: test_solver.o solver_module.o
	$(FC) -o test_solver test_solver.o solver_module.o $(FFLAGS) -llapack -lblas


# Clean the build directory
clean:
	@ rm -f *.o test_solver
