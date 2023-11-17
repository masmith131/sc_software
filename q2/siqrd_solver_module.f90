module siqrd_solver

    ! use solver_gfortran
    use solver  
    implicit none 

    ! integer, parameter :: s_p = kind(1.0) ! single precision 
    ! integer, parameter :: d_p = kind(1.0d0)  ! double precision 
    ! ! following line to change for different precision 
    ! integer, parameter, public :: wp=d_p ! set select_p to the desired precision (sp or dp)
    
    real(wp) :: beta, mu, gamma, alpha, delta 
    real(wp) :: T ! simulation horizon 
    integer (wip):: N ! N +1:  number of grid points in time interval [0,T]
    real(wp) :: S0, I0

    contains

    ! ===============================================================================================
    ! subroutine: setting_parameters
    ! sets the global varibales that correspond to parameters of the siqrd model 
    ! ===============================================================================================
    subroutine setting_parameters(param, x0, arg1, arg2) 
        real(wp), intent(in) :: param(5), x0(2), arg1
        integer(wip), intent(in) :: arg2 
        ! param = (beta, mu, gamma, alpha, delta)
        beta = param(1)  ! infection rate 
        mu = param(2)    ! the rate at which immune people become again susceptible
        gamma = param(3) ! recovery rate 
        alpha = param(4) ! death rate 
        delta = param(5) ! rate at which infected people get quarantined
        S0 = x0(1)
        I0 = x0(2)
        T = arg1 
        N = arg2 

    end subroutine setting_parameters
    
    ! ===============================================================================================
    ! subroutine: func 
    ! computes the time derivatives of s,i,q,r,d, so the vector ! fx = (s'(t),i'(t),q'(t),r'(t),d'(t))
    ! based on the SIQRD models 
    ! ================================================================================================
    subroutine func(xk, fx)
        real(wp), dimension(5), intent(in) :: xk ! xk = (s(t), i(t), q(t), r(t), d(t)), point at which evalutate the function
        real(wp), dimension(5), intent(out) :: fx ! to store output of time derivates of the 5 functions at xk

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
        real(wp), dimension(5), intent(in) :: xk ! xk = (s(t), i(t), q(t), r(t), d(t)), point at which evalutate the jacobian
        real(wp), dimension(5,5), intent(out) :: dfdx  ! dfdx: 5x5 matrix that will store the jacobian 

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
        real(wp), dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real(wp), dimension(5), intent(out) :: xk1 !xk1 the solution at time step k+1
        real(wp), dimension(5) :: fx ! f(x_k)
        call func(xk, fx)
        xk1 = xk + T/N * fx
    end subroutine 

    ! ===========================================================================
    ! subroutine: heun
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Heun's scheme 
    ! ============================================================================    
    subroutine heun(xk, xk1) 
        real(wp), dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real(wp), dimension(5), intent(out) :: xk1 !xk1 the solution at time step k+1
        real(wp), dimension(5) ::  fx, xkt, fxt! f(x_k), x_k + T/N f(x_k), f(x_k + T/N f(x_k))
        call func(xk, fx)
        xkt = xk + T/N*fx
        call func(xkt, fxt)
        xk1 = xk + T/N * (0.5_wp * fx + 0.5_wp * fxt)     
    end subroutine 

    ! ===========================================================================
    ! subroutine: backward
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Euler's backward scheme 
    ! ============================================================================  
    subroutine backward(xk,xk1)
        real(wp), dimension(5), intent(in) :: xk ! x_k the solution at time step k 
        real(wp), dimension(5), intent(out) :: xk1!xk1 the solution at time step k+1
        real(wp), dimension(5) :: fx(5), be(5), prod(5) ! f(x_k+1^(s+1)), corresponds to the last term of newton or the backward error 
        real(wp), dimension(5,5) :: dfdx, Id !dfdx will contain the matrix term and Id the identity 5 by 5 
        real(wp), parameter :: rtol = 1e-15_wp

        integer, parameter :: max_it = 200 !maximum number of iterations of newton for one step
        integer i

        Id = 0.0_wp
        !Identity Matrix
        do i = 1,5
            Id(i,i) = 1.0_wp
        enddo

        ! initial guess set to xk
        xk1 = xk 

        ! iterations of newton's method 
        do i = 1,max_it 
            if(i == max_it) then 
                print *, "Warning the newton method reached its maximum number of iterations !"
                stop 
            endif

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

            !backward error stopping criteria
            if(norm2(be)/norm2(xk1) < rtol) return 

        enddo 
    end subroutine 

    !------------------------------------------------------------------------------------------------------------
    ! f = iq_max(b,s0,i0)
    ! computes the maximum value of Ik + Qk for a given beta (b)
    !-------------------------------------------------------------------------------------------------------------
    
    function iq_max(b, meth) result(f)
        real(wp), intent(in) :: b  ! beta 
        real(wp) :: f              ! will contain the result F(beta)
        character :: meth             ! method used to find solution       
        real(wp), dimension(N+1) :: solk  ! array to store solution at time step k 
        real(wp), dimension(N+1) :: solk1 ! array to store solution at time step k+1
        integer(wip) i

        beta = b        ! computes I+Q as a function of the input b
        
        ! fill in solution at time step 0
        solk  = 0.0_wp
        solk(1) = S0
        solk(2) = I0
        f = I0

        ! get solution at time k+1 from solution at time k using method 
        do i = 2, N+1
            if(meth == 'f') then 
                call forward(solk, solk1)
            elseif(meth == 'h') then 
                call heun(solk, solk1)
            elseif(meth == 'b') then 
                call backward(solk, solk1)
            else 
                print *, "method m is not recognized"
                exit
            endif
            if (f < (solk1(2)+solk1(3))) f =  (solk1(2)+solk1(3))

            solk = solk1
            solk1 = 0.0_wp
        enddo

    end function 


end module siqrd_solver