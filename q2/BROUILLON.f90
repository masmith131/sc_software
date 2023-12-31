
    subroutine backward(param, xk, T, N, xk1)
        real param(5), xk(5), T, xk1(5)
        integer i,N, max_it
        max_it = 200

        ! initial guess set to xk
        xk1 = xk
        do i = 1,max_it 
            call newton(param,xk1,T,N,xk)
            call func(param,xs, fx )
            if (norm2(xk + T/N*fx - xs) < 1.0) exit 
        enddo 

    end subroutine 


        ! subroutine: newton 
    ! computes one step of the newton method 
    subroutine newton(param, xk1, T, N,xk)
        ! xs corresponds to x_(k+1)^(s) 
        ! xs1 corresponds to x_(k+1)^(s+1)
        ! inv corresponds to the second term of the newton method 
        real T 
        real, dimension(5) :: param, xk1, xk, fx, incr
        real, dimension(5,5) :: dfdx, Id
        integer N, i
        
        !Identity Matrix
        do i = 1,5
            Id(i,i) = 1.0
        enddo

        !computing the second term (before inversion)
        call jacob(param,xs, dfdx)  ! compute the jacobian 
        dfdx = dfdx * T/N           ! scale the jacobian 
        dfdx = dfdx - Id            ! substract it the identity

        ! computing the third term and putting it in xk
        call func(param,xk1, fx)
        incr = xk + T/N * fx - xs

        ! solving (dfdx)x = incr to obtain x = (dfdx)^-1 * incr
        call solve(dfdx,incr)
        
        ! final computation 
        xk1 = xk1 - incr

    end subroutine 

    subroutine backward(xk, xk1)
        real, dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real, dimension(5), intent(out) :: xk1 !xk1 the solution at time step k+1
        real, dimension(5) :: fx, res, deltax
        real, dimension(5,5) :: J, I
        integer :: iter, max_iter
        real :: tol, error
        
        max_iter = 100
        tol = 1.0e-6
        I = 0.0
        I(1,1) = 1.0
        I(2,2) = 1.0
        I(3,3) = 1.0
        I(4,4) = 1.0
        I(5,5) = 1.0
        
        ! Initial guess for xk1
        xk1 = xk
        
        do iter = 1, max_iter
            call func(xk1, fx)
            call jacob(xk1, J)
            
            ! Compute the residual
            res = xk + T/N * fx - xk1
            
            ! Newton's step
            ! solve for deltax: (I - T/N * J) * deltax = res
            J = I - T/N * J
            call solve(J, res)
            deltax = res
            
            ! Update the estimate
            xk1 = xk1 + deltax
            
            ! Check convergence
            error = sum(abs(deltax))
            if (error < tol) exit
        end do
    end subroutine
    
    fx(1) = -beta *(xk(2)/(xk(1)+ xk(2)+ xk(4)))* xk(1)+ mu*xk(4)
    fx(2) = (beta* xk(1)/(xk(1) + xk(2) + xk(4)) - gamma - delta - alpha) * xk(2)
    fx(3) = delta* xk(2) - (gamma + alpha)*xk(3)
    fx(4) = gamma* (xk(2) + xk(3)) - mu *xk(4)
    fx(5) = alpha* (xk(2) + xk(3))

    1.00000E+02   5.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  
    9.36103E+01   3.55929E+00  -2.40641E-13   7.11857E+00   7.11857E-01  
    8.94555E+01   2.41064E+00   5.55250E-13   1.19399E+01   1.19399E+00  
    8.67996E+01   1.58330E+00  -2.99018E-13   1.51065E+01   1.51065E+00  
    8.51179E+01   1.02031E+00  -2.51191E-13   1.71471E+01   1.71471E+00  
    8.40589E+01   6.49776E-01  -1.34115E-12   1.84466E+01   1.84466E+00  
    8.33942E+01   4.10771E-01  -3.79546E-13   1.92682E+01   1.92682E+00  
    8.29778E+01   2.58490E-01  -1.18608E-13   1.97852E+01   1.97852E+00  
    8.27173E+01   1.62198E-01  -5.68434E-14   2.01096E+01   2.01096E+00  
    8.25544E+01   1.01595E-01  -1.77636E-14   2.03127E+01   2.03127E+00  
    8.24526E+01   6.35647E-02  -1.54599E-14   2.04399E+01   2.04399E+00  
    8.23890E+01   3.97427E-02   1.26635E-16   2.05194E+01   2.05194E+00  
    8.23492E+01   2.48376E-02   7.47839E-15   2.05690E+01   2.05690E+00  
    8.23244E+01   1.55183E-02   2.33581E-15   2.06001E+01   2.06001E+00  
    8.23089E+01   9.69404E-03  -5.10009E-16   2.06195E+01   2.06195E+00  
    8.22992E+01   6.05508E-03   8.20677E-13   2.06316E+01   2.06316E+00  
    8.22932E+01   3.78187E-03   2.56771E-13   2.06391E+01   2.06391E+00  
    8.22894E+01   2.36197E-03   2.04725E-13   2.06439E+01   2.06439E+00  
    8.22870E+01   1.47514E-03   1.13243E-13   2.06468E+01   2.06468E+00  
    8.22856E+01   9.21261E-04   5.41789E-14   2.06486E+01   2.06486E+00  
    8.22846E+01   5.75345E-04   2.44249E-14   2.06498E+01   2.06498E+00  
    8.22841E+01   3.59312E-04   7.59375E-15   2.06505E+01   2.06505E+00  
    8.22837E+01   2.24395E-04   2.35367E-15   2.06510E+01   2.06510E+00  
    8.22835E+01   1.40137E-04   1.22125E-15   2.06512E+01   2.06512E+00  
    8.22834E+01   8.75172E-05   5.13478E-16   2.06514E+01   2.06514E+00  
    8.22833E+01   5.46554E-05   1.60462E-16   2.06515E+01   2.06515E+00  
    8.22832E+01   3.41328E-05   5.01444E-17   2.06516E+01   2.06516E+00  
    8.22832E+01   2.13163E-05   1.56701E-17   2.06516E+01   2.06516E+00  
    8.22832E+01   1.33122E-05   6.93889E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   8.31359E-06   1.73472E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   5.19191E-06  -6.33693E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   3.24240E-06   8.67362E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   2.02491E-06   0.00000E+00   2.06517E+01   2.06517E+00  
    8.22831E+01   1.26457E-06   0.00000E+00   2.06517E+01   2.06517E+00  
    8.22831E+01   7.89736E-07   0.00000E+00   2.06517E+01   2.06517E+00  
    8.22831E+01   4.93197E-07   0.00000E+00   2.06517E+01   2.06517E+00  
    8.22831E+01   3.08006E-07   3.78419E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   1.92353E-07  -2.46572E-13   2.06517E+01   2.06517E+00  
    8.22831E+01   1.20126E-07  -2.31040E-13   2.06517E+01   2.06517E+00  
    8.22831E+01   7.50198E-08  -1.68366E-13   2.06517E+01   2.06517E+00  
    8.22831E+01   4.68505E-08  -1.12671E-13   2.06517E+01   2.06517E+00  
    8.22831E+01   2.92586E-08  -7.27154E-14   2.06517E+01   2.06517E+00  
    8.22831E+01   1.82722E-08  -4.61462E-14   2.06517E+01   2.06517E+00  
    8.22831E+01   1.14112E-08  -2.90484E-14   2.06517E+01   2.06517E+00  
    8.22831E+01   7.12638E-09  -1.82127E-14   2.06517E+01   2.06517E+00  
    8.22831E+01   4.45049E-09  -1.13964E-14   2.06517E+01   2.06517E+00  
    8.22831E+01   2.77937E-09  -7.12419E-15   2.06517E+01   2.06517E+00  
    8.22831E+01   1.73574E-09  -4.45131E-15   2.06517E+01   2.06517E+00  
    8.22831E+01   1.08399E-09  -2.78057E-15   2.06517E+01   2.06517E+00  
    8.22831E+01   6.76959E-10  -1.73670E-15   2.06517E+01   2.06517E+00  
    8.22831E+01   4.22767E-10  -1.08465E-15   2.06517E+01   2.06517E+00  
    8.22831E+01   2.64022E-10  -6.77397E-16   2.06517E+01   2.06517E+00  
    8.22831E+01   1.64884E-10  -4.23047E-16   2.06517E+01   2.06517E+00  
    8.22831E+01   1.02972E-10  -2.64199E-16   2.06517E+01   2.06517E+00  
    8.22831E+01   6.43067E-11  -1.64995E-16   2.06517E+01   2.06517E+00  
    8.22831E+01   4.01601E-11  -1.03041E-16   2.06517E+01   2.06517E+00  
    8.22831E+01   2.50803E-11  -6.43502E-17   2.06517E+01   2.06517E+00  
    8.22831E+01   1.56629E-11  -4.01873E-17   2.06517E+01   2.06517E+00  
    8.22831E+01   9.78162E-12  -2.50973E-17   2.06517E+01   2.06517E+00  
    8.22831E+01   6.10871E-12  -1.56735E-17   2.06517E+01   2.06517E+00  
    8.22831E+01   3.81495E-12  -9.78824E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   2.38247E-12  -6.11285E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   1.48787E-12  -3.81753E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   9.29190E-13  -2.38408E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   5.80287E-13  -1.48888E-18   2.06517E+01   2.06517E+00  
    8.22831E+01   3.62395E-13  -9.29819E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   2.26319E-13  -5.80680E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   1.41338E-13  -3.62640E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   8.82669E-14  -2.26472E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   5.51235E-14  -1.41434E-19   2.06517E+01   2.06517E+00  
    8.22831E+01   3.44251E-14  -8.83266E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   2.14988E-14  -5.51608E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   1.34262E-14  -3.44484E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   8.38477E-15  -2.15133E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   5.23637E-15  -1.34353E-20   2.06517E+01   2.06517E+00  
    8.22831E+01   3.27016E-15  -8.39045E-21   2.06517E+01   2.06517E+00  
    8.22831E+01   2.04224E-15  -5.23991E-21   2.06517E+01   2.06517E+00  
    8.22831E+01   1.27540E-15  -3.27237E-21   2.06517E+01   2.06517E+00  
    8.22831E+01   7.96498E-16  -2.04363E-21   2.06517E+01   2.06517E+00  
    8.22831E+01   4.97420E-16  -1.27626E-21   2.06517E+01   2.06517E+00  
    8.22831E+01   3.10643E-16  -7.97037E-22   2.06517E+01   2.06517E+00  
    8.22831E+01   1.94000E-16  -4.97757E-22   2.06517E+01   2.06517E+00  
    8.22831E+01   1.21155E-16  -3.10854E-22   2.06517E+01   2.06517E+00  
    8.22831E+01   7.56621E-17  -1.94131E-22   2.06517E+01   2.06517E+00  
    8.22831E+01   4.72516E-17  -1.21237E-22   2.06517E+01   2.06517E+00  
    8.22831E+01   2.95091E-17  -7.57133E-23   2.06517E+01   2.06517E+00  
    8.22831E+01   1.84287E-17  -4.72836E-23   2.06517E+01   2.06517E+00  
    8.22831E+01   1.15089E-17  -2.95291E-23   2.06517E+01   2.06517E+00  
    8.22831E+01   7.18740E-18  -1.84412E-23   2.06517E+01   2.06517E+00  
    8.22831E+01   4.48860E-18  -1.15167E-23   2.06517E+01   2.06517E+00  
    8.22831E+01   2.80317E-18  -7.19226E-24   2.06517E+01   2.06517E+00  
    8.22831E+01   1.75060E-18  -4.49163E-24   2.06517E+01   2.06517E+00  
    8.22831E+01   1.09327E-18  -2.80507E-24   2.06517E+01   2.06517E+00  
    8.22831E+01   6.82756E-19  -1.75179E-24   2.06517E+01   2.06517E+00  
    8.22831E+01   4.26387E-19  -1.09401E-24   2.06517E+01   2.06517E+00  
    8.22831E+01   2.66283E-19  -6.83218E-25   2.06517E+01   2.06517E+00  
    8.22831E+01   1.66296E-19  -4.26676E-25   2.06517E+01   2.06517E+00  
    8.22831E+01   1.03853E-19  -2.66463E-25   2.06517E+01   2.06517E+00  
    8.22831E+01   6.48573E-20  -1.66408E-25   2.06517E+01   2.06517E+00  
    8.22831E+01   4.05040E-20  -1.03924E-25   2.06517E+01   2.06517E+00  
    8.22831E+01   2.52951E-20  -6.49012E-26   2.06517E+01   2.06517E+00  
    8.22831E+01   1.57970E-20  -4.05314E-26   2.06517E+01   2.06517E+00  
    8.22831E+01   9.86537E-21  -2.53122E-26   2.06517E+01   2.06517E+00  
    8.22831E+01   6.16101E-21  -1.58077E-26   2.06517E+01   2.06517E+00  
    8.22831E+01   3.84761E-21  -9.87205E-27   2.06517E+01   2.06517E+00  
    8.22831E+01   2.40287E-21  -6.16518E-27   2.06517E+01   2.06517E+00  
    8.22831E+01   1.50061E-21  -3.85021E-27   2.06517E+01   2.06517E+00  
    8.22831E+01   9.37145E-22  -2.40449E-27   2.06517E+01   2.06517E+00  
    8.22831E+01   5.85256E-22  -1.50163E-27   2.06517E+01   2.06517E+00  
    8.22831E+01   3.65497E-22  -9.37780E-28   2.06517E+01   2.06517E+00  
    8.22831E+01   2.28256E-22  -5.85652E-28   2.06517E+01   2.06517E+00  
    8.22831E+01   1.42548E-22  -3.65745E-28   2.06517E+01   2.06517E+00  
    8.22831E+01   8.90226E-23  -2.28411E-28   2.06517E+01   2.06517E+00  
    8.22831E+01   5.55954E-23  -1.42645E-28   2.06517E+01   2.06517E+00  
    8.22831E+01   3.47199E-23  -8.90829E-29   2.06517E+01   2.06517E+00  
    8.22831E+01   2.16829E-23  -5.56331E-29   2.06517E+01   2.06517E+00  
    8.22831E+01   1.35411E-23  -3.47434E-29   2.06517E+01   2.06517E+00  
    8.22831E+01   8.45656E-24  -2.16975E-29   2.06517E+01   2.06517E+00  
    8.22831E+01   5.28120E-24  -1.35503E-29   2.06517E+01   2.06517E+00  
    8.22831E+01   3.29816E-24  -8.46229E-30   2.06517E+01   2.06517E+00  
    8.22831E+01   2.05973E-24  -5.28478E-30   2.06517E+01   2.06517E+00  
    8.22831E+01   1.28632E-24  -3.30039E-30   2.06517E+01   2.06517E+00  
    8.22831E+01   8.03318E-25  -2.06112E-30   2.06517E+01   2.06517E+00  
    8.22831E+01   5.01679E-25  -1.28719E-30   2.06517E+01   2.06517E+00  
    8.22831E+01   3.13303E-25  -8.03862E-31   2.06517E+01   2.06517E+00  
    8.22831E+01   1.95661E-25  -5.02019E-31   2.06517E+01   2.06517E+00  
    8.22831E+01   1.22192E-25  -3.13515E-31   2.06517E+01   2.06517E+00  
    8.22831E+01   7.63099E-26  -1.95793E-31   2.06517E+01   2.06517E+00  
    8.22831E+01   4.76562E-26  -1.22275E-31   2.06517E+01   2.06517E+00  
    8.22831E+01   2.97617E-26  -7.63616E-32   2.06517E+01   2.06517E+00  
    8.22831E+01   1.85865E-26  -4.76885E-32   2.06517E+01   2.06517E+00  
    8.22831E+01   1.16074E-26  -2.97819E-32   2.06517E+01   2.06517E+00  
    8.22831E+01   7.24894E-27  -1.85991E-32   2.06517E+01   2.06517E+00  
    8.22831E+01   4.52703E-27  -1.16153E-32   2.06517E+01   2.06517E+00  
    8.22831E+01   2.82717E-27  -7.25384E-33   2.06517E+01   2.06517E+00  
    8.22831E+01   1.76559E-27  -4.53009E-33   2.06517E+01   2.06517E+00  
    8.22831E+01   1.10263E-27  -2.82908E-33   2.06517E+01   2.06517E+00  
    8.22831E+01   6.88601E-28  -1.76679E-33   2.06517E+01   2.06517E+00  
    8.22831E+01   4.30038E-28  -1.10337E-33   2.06517E+01   2.06517E+00  
    8.22831E+01   2.68562E-28  -6.89067E-34   2.06517E+01   2.06517E+00  
    8.22831E+01   1.67720E-28  -4.30329E-34   2.06517E+01   2.06517E+00  
    8.22831E+01   1.04742E-28  -2.68744E-34   2.06517E+01   2.06517E+00  
    8.22831E+01   6.54126E-29  -1.67833E-34   2.06517E+01   2.06517E+00  
    8.22831E+01   4.08507E-29  -1.04813E-34   2.06517E+01   2.06517E+00  
    8.22831E+01   2.55117E-29  -6.54569E-35   2.06517E+01   2.06517E+00  
    8.22831E+01   1.59323E-29  -4.08784E-35   2.06517E+01   2.06517E+00  
    8.22831E+01   9.94984E-30  -2.55289E-35   2.06517E+01   2.06517E+00  
    8.22831E+01   6.21377E-30  -1.59430E-35   2.06517E+01   2.06517E+00  
    8.22831E+01   3.88055E-30  -9.95658E-36   2.06517E+01   2.06517E+00  
    8.22831E+01   2.42344E-30  -6.21797E-36   2.06517E+01   2.06517E+00  
    8.22831E+01   1.51346E-30  -3.88318E-36   2.06517E+01   2.06517E+00 



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


siqrd_solver_module.o: siqrd_solver_module.f90 solver.mod
	$(FC) -o siqrd_solver_module.o -c siqrd_solver_module.f90

siqrd.o: siqrd.f90 siqrd_solver_module.o
	$(FC) $(FFLAGS) -o siqrd.o -c siqrd.f90

siqrd: siqrd.o siqrd_solver_module.o
	$(FC) -o siqrd siqrd.o solver_module.o siqrd_solver_module.o $(FFLAGS) -llapack -lblas

clean:
	@ rm -f siqrd.o siqrd_solver_module.o solver_module.o test_solver.o siqrd_solver.mod siqrd test_solver



    