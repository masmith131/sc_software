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



CORRECTION 
LEFT TO DO FOR PART 1 
- find out if kind needs to be put for integers if when we use double precision we also need double integer precision 
- make sure everything is defined in correct precision 
- check comments in report about this part (that everything is alligned with what is done in code)

PROGRESS 
thoughts for function iq 
- should method be global var/ argument / defined in the function ? 
- is the beta0 the one given in parameters.in ? or not necessarily 
- should i0 and s0 be global varibales ? 
next up: 
- in tweaking compute f(beta) and f(beta+deltab) and then update the b0 ... 
- think about deltab and about stopping criteria 