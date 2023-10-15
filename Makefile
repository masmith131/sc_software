FC=gfortran

FFLAGS_g95      = -g -pedantic -Wall -fbounds-check -ftrace=full
FFLAGS_gfortran = -g -pedantic -Wall -Wimplicit-interface -Wunderflow -fbounds-check -fimplicit-none
FFLAGS_ifort    = -g -debug full -implicitnone -check -warn -free -Tf
FFLAGS_nagfor   = -g -gline -u -info -colour


# Select the right flags for the current compiler
FFLAGS=$(FFLAGS_$(FC))

all: fortran_start
	./fortran_start

fortran_start.o: fortran_start.f90 
	gfortran -o fortran_start.o -c fortran_start.f90

fortran_start: fortran_start.o 
	gfortran -o fortran_start fortran_start.o

clean:
	@ rm -f *.o *.mod fortran_start