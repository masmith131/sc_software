module siqrd_solver

    use solver_gfortran 
    implicit none 

    real :: beta, mu, gamma, alpha, delta
    real, parameter, public :: T = 30.0 ! simulation horizon 
    integer, parameter, public  :: N = 150 ! N +1:  number of grid points in time interval [0,T]

    contains

    ! ===============================================================================================
    ! subroutine: setting_parameters
    ! sets the global varibales that correspond to parameters of the siqrd model 
    ! ===============================================================================================
    subroutine setting_parameters(param)
        real, intent(in) :: param(5)
        ! param = (beta, mu, gamma, alpha, delta)
        beta = param(1)  ! infection rate 
        mu = param(2)    ! the rate at which immune people become again susceptible
        gamma = param(3) ! recovery rate 
        alpha = param(4) ! rate at which infected people get tested
        delta = param(5) ! death rate 
    end subroutine setting_parameters
    
    ! ===============================================================================================
    ! subroutine: func 
    ! computes the time derivatives of s,i,q,r,d, so the vector ! fx = (s'(t),i'(t),q'(t),r'(t),d'(t))
    ! based on the SIQRD models 
    ! ================================================================================================
    subroutine func(xk, fx)
        real, dimension(5), intent(in) :: xk ! xk = (s(t), i(t), q(t), r(t), d(t)), point at which evalutate the function
        real, dimension(5), intent(out) :: fx ! to store output of time derivates of the 5 functions at xk

        ! xk = (s(t), i(t), q(t), r(t), d(t))
        ! fx = (s'(t),i'(t),q'(t),r'(t),d'(t))
        fx(1) = -beta *(xk(2)/(xk(1)+ xk(2)+ xk(4)))* xk(1)+ mu*xk(4)
        fx(2) = (beta* xk(1)/(xk(1) + xk(2) + xk(4)) - gamma - delta - alpha) * xk(2)
        fx(3) = delta* xk(2) - (gamma + alpha)*xk(3)
        fx(4) = gamma* (xk(2) + xk(3)) - mu *xk(4)
        fx(5) = alpha* (xk(2) + xk(3))

    end subroutine func 

    ! ===============================================================================================
    ! subroutine: jacob
    ! computes the jacobian of function f(t) = (s'(t),i'(t),q'(t),r'(t),d'(t))
    ! ===============================================================================================

    subroutine jacob(xk, dfdx)
        real, dimension(5), intent(in) :: xk ! xk = (s(t), i(t), q(t), r(t), d(t)), point at which evalutate the jacobian
        real , dimension(5,5), intent(out) :: dfdx  ! dfdx: 5x5 matrix that will store the jacobian 

        dfdx = 0.0
        ! computing all non zero entries of the jacobian 
        dfdx(1,1) = - beta*xk(2) *((xk(2)+xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(1,2) = - beta*xk(1) *((xk(1) + xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(1,4) = (beta*xk(1)*xk(2))/((xk(1)+xk(2)+xk(4))**2) + mu
        dfdx(2,1) = beta*xk(2) *((xk(2)+xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(2,2) = beta*xk(1) *((xk(1) + xk(4))/(xk(1)+xk(2)+xk(4))**2) - gamma - delta - alpha
        dfdx(2,4) = (-beta*xk(1)*xk(2))/((xk(1)+xk(2)+xk(4))**2)
        dfdx(3,2) = delta 
        dfdx(3,3) = - gamma - alpha
        dfdx(4,2) = gamma
        dfdx(4,3) = gamma
        dfdx(4,4) = - mu 
        dfdx(5,2) = alpha 
        dfdx(5,3) = alpha
    end subroutine jacob

    ! ===========================================================================
    ! subroutine: forward
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Euler's forward scheme 
    ! ============================================================================
    subroutine forward(xk,xk1)
        real, dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real, dimension(5), intent(out) :: xk1 !xk1 the solution at time step k+1
        real, dimension(5) :: fx ! f(x_k)
        call func(xk, fx)
        xk1 = xk + T/N * fx
    end subroutine 

    ! ===========================================================================
    ! subroutine: heun
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Heun's scheme 
    ! ============================================================================    
    subroutine heun(xk, xk1) 
        real, dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real, dimension(5), intent(out) :: xk1 !xk1 the solution at time step k+1
        real, dimension(5) ::  fx, xkt, fxt! f(x_k), x_k + T/N f(x_k), f(x_k + T/N f(x_k))
        call func(xk, fx)
        xkt = xk + T/N*fx
        call func(xkt, fxt)
        xk1 = xk + T/N * (0.5 * fx + 0.5 * fxt)     
    end subroutine 

    ! ===========================================================================
    ! subroutine: backward
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Euler's backward scheme 
    ! ============================================================================  
    subroutine backward(xk,xk1)
        real, dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real, dimension(5), intent(out) :: xk1!xk1 the solution at time step k+1
        real, dimension(5) :: fx(5), be(5), prod(5) ! f(x_k+1^(s+1)), corresponds to the last term of newton or the backward error 
        real, dimension(5,5) :: dfdx, Id !dfdx will contain the matrix term and Id the identity 5 by 5 

        integer, parameter :: max_it = 100 !maximum number of iterations of newton for one step
        integer i

        !Identity Matrix
        do i = 1,5
            Id(i,i) = 1.0
        enddo

        ! initial guess set to xk
        xk1 = xk 

        ! iterations of newton's method 
        do i = 1,max_it 
            if(i == max_it) print *, "Warning the newton method reached its maximum number of iterations !"

            !computing the second term (before inversion)
            call jacob(xk1, dfdx)  ! compute the jacobian 
            dfdx = dfdx * T/N           ! scale the jacobian 
            dfdx = dfdx - Id            ! substract it the identity

            ! computing the third term 
            call func(xk1, fx)
            be = xk + T/N * fx - xk1
            prod = be 
            
            ! solving (dfdx)x = prod to obtain x = (dfdx)^-1 * prod and result x will be in prod 
            call solve(dfdx,prod)
        
            ! final computation 
            xk1 = xk1 - prod ! now xk1 contains x_(k+1)^(s+1)
            
            print *, xk1 

            !backward error stopping criteria
            if(norm2(be) < 1.0e-6) exit

        enddo 
    end subroutine 


end module siqrd_solver