FC=gfortran

FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -g -gline -u -info -colour

# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))

all: infections_exp1 infections_exp2 infections_exp3

exp1: infections_exp1
	./infections_exp1

exp2: infections_exp2
	./infections_exp2

exp3: infections_exp3
	./infections_exp3

infections_module.o: infections_module.f90
	$(FC) $(FFLAGS) -o infections_module.o -c infections_module.f90

infections_exp1.o: infections_exp1.f90 infections_module.o
	$(FC) $(FFLAGS) -o infections_exp1.o -c infections_exp1.f90

infections_exp2.o: infections_exp2.f90 infections_module.o
	$(FC) $(FFLAGS) -o infections_exp2.o -c infections_exp2.f90

infections_exp3.o: infections_exp3.f90 infections_module.o
	$(FC) $(FFLAGS) -o infections_exp3.o -c infections_exp3.f90

infections_exp1: infections_exp1.o infections_module.o
	$(FC) $(FFLAGS) -o infections_exp1 infections_exp1.o infections_module.o 

infections_exp2: infections_exp2.o infections_module.o
	$(FC) $(FFLAGS) -o infections_exp2 infections_exp2.o infections_module.o 

infections_exp3: infections_exp3.o infections_module.o
	$(FC) $(FFLAGS) -o infections_exp3 infections_exp3.o infections_module.o 


clean:
	@ rm -f infections_exp1.o infections_module.o *.mod infections_exp1
