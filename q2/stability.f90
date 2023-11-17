program stability
    use siqrd_solver

    implicit none

    real(wp), dimension(7) :: input ! stores the input read from 'parameters.in'
    real(wp), dimension(5) :: x     ! vector siqrd 
    real(wp), dimension(5,5) :: J   ! jacobian 
    real(wp), dimension(5) :: lambda !eigenvalues 
    real(wp) :: eps = 1e-15_wp 
    integer i
    logical stable  

    ! getting parameters and s0, i0 from the input file 
    open(unit = 32,file = 'eigenvalues1.in')
    read(32,*) input
    close(unit = 32)
    ! setting parameters for the numerical method method 
    call setting_parameters(input(:5), input(6:), 1.0_wp, 1_wip)  

    x = 0.0_wp 
    x(1) = input(6)
    x(2) = input(7)
    J = 0.0_wp
    call jacob(x, J)
    call get_eigenvalues(J, lambda)
    ! check all eigenvalues are negative
    stable = .True.
    do i = 1,5
        if (lambda(i) > eps) then 
            stable = .False.
            exit
        endif 
    enddo 
    print *, lambda
    if(stable) then 
        print *, "The equilibrium is stable"
    else 
        print *, "The equilibrium is not stable"
    endif
    
end program stability
