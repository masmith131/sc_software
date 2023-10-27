FC=ifort

all: infections_exp1 infections_exp2 

exp1: infections_exp1
	./infections_exp1

exp2: infections_exp2
	./infections_exp2


infections_module.o: infections_module.f90
	$(FC) -o infections_module.o -c infections_module.f90

infections_exp1.o: infections_exp1.f90 infections_module.o
	$(FC) -o infections_exp1.o -c infections_exp1.f90

infections_exp2.o: infections_exp2.f90 infections_module.o
	$(FC) -o infections_exp2.o -c infections_exp2.f90


infections_exp1: infections_exp1.o infections_module.o
	$(FC) -o infections_exp1 infections_exp1.o infections_module.o 

infections_exp2: infections_exp2.o infections_module.o
	$(FC) -o infections_exp2 infections_exp2.o infections_module.o 


clean:
	@ rm -f *.o *.mod infections_exp1 infections_exp2