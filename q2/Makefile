FC=gfortran

FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -g -gline -u -info -colour

# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))

all: siqrd 
	./siqrd

siqrd_solver_module.o: siqrd_solver_module.f90 solver_gfortran.mod
	$(FC) $(FFLAGS) -o siqrd_solver_module.o -c siqrd_solver_module.f90

siqrd.o: siqrd.f90 siqrd_solver_module.o
	$(FC) $(FFLAGS) -o siqrd.o -c siqrd.f90

siqrd: siqrd.o siqrd_solver_module.o
	$(FC) $(FFLAGS) -o siqrd siqrd.o siqrd_solver_module.o solver_gfortran.o

clean:
	@ rm -f siqrd.o siqrd_solver_module.o siqrd_solver.mod siqrd 