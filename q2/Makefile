FC=gfortran

# Flags for g95 that you might also want to use with gfortran
FFLAGS=-g -pedantic -Wall -fbounds-check 

all: siqrd 
	./siqrd 30.0 150

test: test_solver
	./test_solver

solver_module.o: solver_module.f90
	$(FC) -o solver_module.o -c solver_module.f90

test_solver.o: test_solver.f90 solver_module.o
	$(FC) $(FFLAGS) -c test_solver.f90 -o test_solver.o

test_solver: test_solver.o solver_module.o
	$(FC) -o test_solver test_solver.o solver_module.o $(FFLAGS) -llapack -lblas

siqrd_solver_module.o: siqrd_solver_module.f90 solver_module.o
	$(FC) $(FFLAGS) -c siqrd_solver_module.f90 -o siqrd_solver_module.o -llapack -lblas

siqrd.o: siqrd.f90 siqrd_solver_module.o solver_module.o
	$(FC) $(FFLAGS) -c siqrd.f90 -o siqrd.o

siqrd: siqrd.o siqrd_solver_module.o solver_module.o
	$(FC) -o siqrd siqrd.o siqrd_solver_module.o solver_module.o $(FFLAGS) -llapack -lblas


# siqrd_solver_module.o: siqrd_solver_module.f90 solver_gfortran.mod
# 	$(FC) -o siqrd_solver_module.o -c siqrd_solver_module.f90

# siqrd.o: siqrd.f90 siqrd_solver_module.o
# 	$(FC) -o siqrd.o -c siqrd.f90

# siqrd: siqrd.o siqrd_solver_module.o
# 	$(FC) -o siqrd siqrd.o siqrd_solver_module.o solver_gfortran.o


clean:
	@ rm -f siqrd.o siqrd_solver_module.o solver_module.o siqrd_solver.mod siqrd test_solver.o test_solver