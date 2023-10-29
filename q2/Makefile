FC=gfortran

all: siqrd 
	./siqrd

siqrd_solver_module.o: siqrd_solver_module.f90 solver_gfortran.mod
	$(FC) -o siqrd_solver_module.o -c siqrd_solver_module.f90

siqrd.o: siqrd.f90 siqrd_solver_module.o
	$(FC) -o siqrd.o -c siqrd.f90

siqrd: siqrd.o siqrd_solver_module.o
	$(FC) -o siqrd siqrd.o siqrd_solver_module.o solver_gfortran.o

clean:
	@ rm -f siqrd.o siqrd_solver_module.o siqrd_solver.mod siqrd 