
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