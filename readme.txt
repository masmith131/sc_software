The code that answers question 1 is in ./q1
To run it go to ./q1

Commands:
make exp1: will run code for first part of question 1, where the user is prompt
make exp2: will run code for second part of question 1 where parameters are imposed 
make clean: removes executables and object files 

!! to change the rate go to line 16 of 'infections_exp2.f90'
!! to change the precision go to line 8 of 'infections_module'
!! to change the compiler change line 1 of 'Makefile'

The code that answers question 2 is in ./q2
To run it go to ./q2
Commands:

make: will run code and print results to stdout
make clean: removes executables and object files 
!! to change solving method go to line 6 of 'siqrd.f90'
!! to change parameters change them in 'parameters.in' 
!! to change number of steps and time horizon go to line 7 and 8 of 'siqrd_solver_module'
