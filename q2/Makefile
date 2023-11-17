FC=gfortran

# Flags for g95 that you might also want to use with gfortran
FFLAGS=-g -pedantic -Wall -fbounds-check 
OPT_FLAGS=  #-Ofast #-ffast-math -march=native 

all: question2 
	./question2 1000 100.0

tweak: tweaking
	./tweaking 1000 100.0 0.1 20000
valgrind: tweaking
	valgrind --leak-check=yes --show-leak-kinds=all ./tweaking 1000 100.0 0.1 20000


test: test_solver
	./test_solver

solver_module.o: solver_module.f90
	$(FC) $(OPT_FLAGS) -o solver_module.o -c solver_module.f90

test_solver.o: test_solver.f90 solver_module.o
	$(FC) $(FFLAGS) -c test_solver.f90 -o test_solver.o

test_solver: test_solver.o solver_module.o
	$(FC) -o test_solver test_solver.o solver_module.o $(FFLAGS) -llapack -lblas

siqrd_solver_module.o: siqrd_solver_module.f90 solver_module.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) -c siqrd_solver_module.f90 -o siqrd_solver_module.o 

question2.o: question2.f90 siqrd_solver_module.o solver_module.o
	$(FC) $(FFLAGS) -c question2.f90 -o question2.o

question2: question2.o siqrd_solver_module.o solver_module.o
	$(FC) -o question2 question2.o siqrd_solver_module.o solver_module.o $(FFLAGS) -llapack -lblas

tweaking.o: tweaking.f90 siqrd_solver_module.o solver_module.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) -c tweaking.f90 -o tweaking.o

tweaking: tweaking.o siqrd_solver_module.o solver_module.o
	$(FC) -o tweaking tweaking.o siqrd_solver_module.o solver_module.o $(FFLAGS) -llapack -lblas


# siqrd_solver_module.o: siqrd_solver_module.f90 solver_gfortran.mod
# 	$(FC) -o siqrd_solver_module.o -c siqrd_solver_module.f90

# siqrd.o: siqrd.f90 siqrd_solver_module.o
# 	$(FC) -o siqrd.o -c siqrd.f90

# siqrd: siqrd.o siqrd_solver_module.o
# 	$(FC) -o siqrd siqrd.o siqrd_solver_module.o solver_gfortran.o

clean:
	@ rm -f tweaking.o question2.o siqrd_solver_module.o solver_module.o siqrd_solver.mod question2 test_solver.o test_solver tweaking