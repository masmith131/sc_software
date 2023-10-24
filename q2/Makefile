FC=gfortran

FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -g -gline -u -info -colour

# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))

all: siqrd 
	./siqrd

models_module.o: models_module.f90
	$(FC) $(FFLAGS) -o models_module.o -c models_module.f90

siqrd.o: siqrd.f90 models_module.o
	$(FC) $(FFLAGS) -o siqrd.o -c siqrd.f90

siqrd: siqrd.o models_module.o
	$(FC) $(FFLAGS) -o siqrd siqrd.o models_module.o 

clean:
	@ rm -f *.o *.mod siqrd 