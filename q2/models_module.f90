module models 

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
        dfdx(1,4) = (beta*xk(1)*xk(2))/((xk(1)+xk(2)+xk(4))**2) + mu
        dfdx(2,1) = beta*xk(2) *((xk(2)+xk(4))/(xk(1)+xk(2)+xk(4))**2)
        dfdx(2,2) = beta*xk(1) *((xk(1) + xk(4))/(xk(1)+xk(2)+xk(4))**2) - gamma - delta - alpha
        dfdx(2,4) = (beta*xk(1)*xk(2))/((xk(1)+xk(2)+xk(4))**2)
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
    subroutine newton(param, xs, T, N, xs1, xk, fx) 
        ! xs corresponds to x_(k+1)^(s) 
        ! xs1 corresponds to x_(k+1)^(s)
        ! inv corresponds to the second term of the newton method 
        real param(5), xs(5), xs1(5), T, xk(5),fx(5), dfdx(5,5), inv(5,5)
        integer N, i
        
        !computing the second term 
        call jacob(param,xs, dfdx)  ! compute the jacobian 
        dfdx = dfdx * T/N           ! scale the jacobian 
        do i = 1,5                  ! substract it the identity
            dfdx(i,i) = dfdx(i,i) -1
        enddo

        ! inverting 
        do i = 1,5
            inv(:,i) = 0.0
            inv(i,i) = 1.0
            call solve(dfdx, inv(:,i))
        enddo

        ! computing the third term 
        do i = 1,5
            xk(i) = xk(i) + T/N *fx(i) - xs(i)
        enddo

        ! multiplying the second and third term 

        ! final computation 
        xs1 = 0.0


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

    subroutine backward 
    end subroutine 


end module models