module siqrd_solver

    use solver_gfortran 
    implicit none 
    contains
    
    ! subroutine: func 
    ! computes the time derivatives of s,i,q,r,d
    ! based on the SIQRD models 
    ! func computes vector fx based on the 5 parameters and 5 functions of time
    subroutine func(param, xk, fx)
        real beta, mu, gamma, alpha, delta, param(5), xk(5), fx(5)
        ! xk = (s(t), i(t), q(t), r(t), d(t))
        ! param = (beta, mu, gamma, alpha, delta)
        beta = param(1)
        mu = param(2)
        gamma = param(3)
        alpha = param(4)
        delta = param(5)

        ! fx = (s'(t),i'(t),q'(t),r'(t),d'(t))
        fx(1) = -beta *(xk(2)/(xk(1)+ xk(2)+ xk(4)))* xk(1)+ mu*xk(4)
        fx(2) = (beta* xk(1)/(xk(1) + xk(2) + xk(4)) - gamma - delta - alpha) * xk(2)
        fx(3) = delta* xk(2) - (gamma + alpha)*xk(3)
        fx(4) = gamma* (xk(2) + xk(3)) - mu *xk(4)
        fx(5) = alpha* (xk(2) + xk(3))
    end subroutine func 


    ! subroutine: jacob
    ! computes the jacobian of function f(t) = (s'(t),i'(t),q'(t),r'(t),d'(t))
    subroutine jacob(param, xk, dfdx)
        ! param = (beta, mu, gamma, alpha, delta), parameters of siqrd model
        ! xk = (s(t), i(t), q(t), r(t), d(t)), point at which we want to evalutate the jacobian
        ! dfdx: 5x5 matrix that will stock the jacobian 
        real beta, mu, gamma, alpha, delta, param(5), xk(5), dfdx(5,5)
        beta = param(1)
        mu = param(2)
        gamma = param(3)
        alpha = param(4)
        delta = param(5)

        dfdx = 0.0
        ! computing all non zero entries of the jacobian 
        dfdx(1,1) = - beta*xk(2) *((xk(2)+xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(1,2) = - beta*xk(1) *((xk(1) + xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(1,4) = (-beta*xk(1)*xk(2))/((xk(1)+xk(2)+xk(4))**2) + mu
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


    ! subroutine: newton 
    ! computes one step of the newton method 
    subroutine newton(param, xk1, T, N,xk)
        ! xk corresponds to x_k
        ! xk1 corresponds to x_(k+1)^(s)
        real T 
        real, dimension(5) :: param, xk1, xk, fx, incr
        real, dimension(5,5) :: dfdx, Id
        integer N, i
        
        !Identity Matrix
        do i = 1,5
            Id(i,i) = 1.0
        enddo

        !computing the second term (before inversion)
        call jacob(param,xk1, dfdx)  ! compute the jacobian 
        dfdx = dfdx * T/N           ! scale the jacobian 
        dfdx = dfdx - Id            ! substract it the identity

        ! computing the third term and putting it in xk
        call func(param,xk1, fx)
        incr = xk + T/N * fx - xk1

        ! solving (dfdx)x = incr to obtain x = (dfdx)^-1 * incr
        call solve(dfdx,incr)
        
        ! final computation 
        xk1 = xk1 - incr

    end subroutine 

    ! subroutine: forward
    ! computes values of vector x at time k+1 thanks to values at time k
    ! following Euler's forward scheme 
    subroutine forward(param, xk, T, N, xk1)
        real param(5), xk(5), T, xk1(5), fx(5)
        integer i,N
        call func(param, xk, fx)
        do i = 1,5
            xk1(i) = xk(i) + T/N * fx(i)            
        enddo
    end subroutine 


    subroutine heun(param, xk, T, N, xk1) 
        real param(5), xk(5), T, xk1(5), fx(5), xkt(5), fxt(5)
        integer i,N
        call func(param, xk, fx)
        do i = 1,5
            xkt(i) = xk(i) + T/N*fx(i)
        enddo

        call func(param, xkt, fxt)
        do i = 1,5
            xk1(i) = xk(i) + T/N * (0.5 * fx(i) + 0.5 * fxt(i))     
        enddo
    end subroutine 

    subroutine backward(param, xk, T, N, xk1)
        real param(5), xk(5), T, xk1(5), fxk1(5), err
        integer i,N, max_it
        max_it = 200

        ! initial guess set to xk
        xk1 = (/80.0, 0.0, 0.0, 20.0, 5.0/)
        do i = 1,max_it 
            call newton(param,xk1,T,N,xk)
            ! computing backward error 
            call func(param, xk1, fxk1)
            err = norm2(xk - T/N * fxk1 - xk1)
            if (err < 0.1) exit 
        enddo 

    end subroutine 


end module siqrd_solver